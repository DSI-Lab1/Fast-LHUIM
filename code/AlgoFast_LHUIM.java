package algo;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * 
 * @author RenGuibin
 */
public class AlgoFast_LHUIM {

	/** the set of high-utility itemsets */
	private Itemsets highUtilityItemsets;

	/** object to write the output file */
	BufferedWriter writer = null;

	/** the number of high-utility itemsets found (for statistics) */
	private long patternCount;

	/** the start time and end time of the last algorithm execution */
	long startTimestamp;
	long endTimestamp;

	/** the minutil threshold */
	int lminUtil;

	/**
	 * if this variable is set to true, some debugging information will be shown
	 */
	final boolean DEBUG = false;
	
	/**
	 * If true, use optimization 1, otherwise use optimization 2
	 */
	final boolean OPT = true;

	/**
	 * The following variables are the utility-bins array // Recall that each
	 * bucket correspond to an item
	 */
	/** utility bin array for sub-tree utility */
	private SU[] utilityBinArraySU;
	/** utility bin array for local utility */
	private LU[] utilityBinArrayLU;
	/** a temporary buffer */
	private int[] temp = new int[500];

	/** The total time spent for performing intersections */
	long timeIntersections;
	/** The total time spent for performing database reduction */
	long timeDatabaseReduction;
	/** The total time spent for identifying promising items */
	long timeIdentifyPromisingItems;
	long windowConsumer,windowConsumerLU;
	/** The total time spent for sorting */
	long timeSort;
	/** The total time spent for binary search */
	long timeBinarySearch;

	/** an array that map an old item name to its new name */
	public int[] oldNameToNewNames;
	/** an array that map a new item name to its old name */
	public int[] newNamesToOldNames;
	/** the number of new items */
	int newItemCount;

	/** if true, transaction merging will be performed by the algorithm */
	boolean activateTransactionMerging = false;

	/** A parameter for transaction merging */
	final int MAXIMUM_SIZE_MERGING = 1000;

	/** number of times a transaction was read */
	long transactionReadingCount;
	/** number of merges */
	public long mergeCount;
	private int minLen;
	/** number of itemsets from the search tree that were considered */
	private long candidateCount;
	/** If true, sub-tree utility pruning will be performed */
	private boolean activateSubtreeUtilityPruning;

	/**
	 * Constructor
	 */
	public AlgoFast_LHUIM() {

	}

	/**
	 * Run the algorithm
	 * 
	 * @param lminUtil
	 *            the minimum utility threshold (a positive integer)
	 * @param inputPath
	 *            the input file path
	 * @param outputPath
	 *            the output file path to save the result or null if to be kept
	 *            in memory
	 * @param activateTransactionMerging
	 * @param activateSubtreeUtilityPruning
	 * @param maximumTransactionCount
	 * @return the itemsets or null if the user choose to save to file
	 * @throws IOException
	 *             if exception while reading/writing to file
	 */
	public Itemsets runAlgorithm(int lminUtil, int minLen, String inputPath,
			String outputPath, boolean activateTransactionMerging,
			int maximumTransactionCount, boolean activateSubtreeUtilityPruning)
			throws IOException {

		// reset variables for statistics
		mergeCount = 0;
		transactionReadingCount = 0;
		timeIntersections = 0;
		timeDatabaseReduction = 0;

		// save parameters about activating or not the optimizations
		this.activateTransactionMerging = activateTransactionMerging;
		this.activateSubtreeUtilityPruning = activateSubtreeUtilityPruning;

		// record the start time
		startTimestamp = System.currentTimeMillis();

		// read the input file
		Dataset dataset = new Dataset(inputPath, maximumTransactionCount);

		// save minUtil value selected by the user
		this.lminUtil = lminUtil;
		this.minLen = minLen;
		// if the user choose to save to file
		// create object for writing the output file
		if (outputPath != null) {
			writer = new BufferedWriter(new FileWriter(outputPath));
		} else {
			// if the user choose to save to memory
			writer = null;
			this.highUtilityItemsets = new Itemsets("Itemsets");
		}

		// reset the number of itemset found
		patternCount = 0;

		// reset the memory usage checking utility
		MemoryLogger.getInstance().reset();
		// if in debug mode, show the initial database in the console
		if (DEBUG) {
			System.out.println("===== Initial database === ");
			System.out.println(dataset.toString());
		}

		// Scan the database using utility-bin array to calculate the TWU
		// of each item
		useUtilityBinArrayToCalculateLocalUtilityFirstTime(dataset);

		// if in debug mode, show the TWU calculated using the utility-bin array

		// Now, we keep only the promising items (those having a twu >= minutil)
		List<Integer> itemsToKeep = new ArrayList<Integer>();

		for (int j = 0; j < utilityBinArrayLU.length; j++) {
			if (utilityBinArrayLU[j] != null
					&& utilityBinArrayLU[j].sumlu >= lminUtil) {
				generatePeriodOfLU(utilityBinArrayLU[j], lminUtil, minLen);
				if (utilityBinArrayLU[j].luPeriod.size() > 0) {
					itemsToKeep.add(j);
				}
			}

		}
		if (DEBUG) {
			System.out.println("===== LU OF SINGLE ITEMS === ");
			for (int i = 1; i < utilityBinArrayLU.length; i++) {
				System.out.println("item : " + i);

				for (PairLU pair : utilityBinArrayLU[i].pairLU) {
					System.out.println("ts : " + pair.timestamp + " lu : "
							+ pair.luOftimestamp);
				}
				if (utilityBinArrayLU[i] != null)
					System.out.println("luPeriod : "
							+ utilityBinArrayLU[i].luPeriod
							+ utilityBinArrayLU[i].luPeriod.size());
			}
		}

		// Sort promising items according to the increasing order of TWU
		insertionSort(itemsToKeep, utilityBinArrayLU);
		if (DEBUG) {
			System.out.println("---------------itemsToKeep-------------");
			System.out.println(itemsToKeep.toString());
			System.out
					.println("---------------ts lu luPeriod of single item-------------------");
			for (int i = 1; i < utilityBinArrayLU.length; i++) {
				System.out.println("item : " + i + " sumlu: "
						+ utilityBinArrayLU[i].sumlu);
				System.out.print("luPeriod : ");
				for (Period period : utilityBinArrayLU[i].luPeriod) {
					System.out.print("[" + period.start + "," + period.end
							+ "]");
				}
				System.out.println();
				for (PairLU pair : utilityBinArrayLU[i].pairLU) {
					System.out.println("ts : " + pair.timestamp + " lu : "
							+ pair.luOftimestamp);
				}

			}
		}
		// Rename promising items according to the increasing order of TWU.
		// This will allow very fast comparison between items later by the
		// algorithm
		// This structure will store the new name corresponding to each old name
		oldNameToNewNames = new int[dataset.getMaxItem() + 1];
		// This structure will store the old name corresponding to each new name
		newNamesToOldNames = new int[dataset.getMaxItem() + 1];
		// We will now give the new names starting from the name "1"
		int currentName = 1;
		// For each item in increasing order of TWU
		for (int j = 0; j < itemsToKeep.size(); j++) {
			// get the item old name
			int item = itemsToKeep.get(j);
			// give it the new name
			oldNameToNewNames[item] = currentName;
			// remember its old name
			newNamesToOldNames[currentName] = item;
			// replace its old name by the new name in the list of promising
			// items
			itemsToKeep.set(j, currentName);
			// increment by one the current name so that
			currentName++;
		}

		// remember the number of promising item
		newItemCount = itemsToKeep.size();
		// initialize the utility-bin array for counting the subtree utility
		utilityBinArraySU = new SU[newItemCount + 1];

		// if in debug mode, print to the old names and new names to the console
		// to check if they are correct
		if (DEBUG) {
			System.out.println(itemsToKeep);
			System.out.println(Arrays.toString(oldNameToNewNames));
			System.out.println(Arrays.toString(newNamesToOldNames));
		}

		// We now loop over each transaction from the dataset
		// to remove unpromising items
		for (int i = 0; i < dataset.getTransactions().size(); i++) {
			// Get the transaction
			Transaction transaction = dataset.getTransactions().get(i);

			// Remove unpromising items from the transaction and at the same
			// time
			// rename the items in the transaction according to their new names
			// and sort the transaction by increasing TWU order
			transaction.removeUnpromisingItems(oldNameToNewNames);
		}

		// Now we will sort transactions in the database according to the
		// proposed
		// total order on transaction (the lexicographical order when
		// transactions
		// are read backward).
		long timeStartSorting = System.currentTimeMillis();

		// record the total time spent for sorting
		timeSort = System.currentTimeMillis() - timeStartSorting;

		// if in debug mode, print the database after sorting and removing
		// promising items
		if (DEBUG) {
			System.out
					.println("===== Database without unpromising items and sorted by TWU increasing order === ");
			System.out.println(dataset.toString());
		}

		// Use an utility-bin array to calculate the sub-tree utility of each
		// item
		useUtilityBinArrayToCalculateSubtreeUtilityFirstTime(dataset);

		// Calculate the set of items that pass the sub-tree utility pruning
		// condition
		List<Integer> itemsToExplore = new ArrayList<Integer>();
		// if subtree utility pruning is activated
		if (activateSubtreeUtilityPruning) {
			// for each item
			for (int i = 0; i < utilityBinArraySU.length; i++) {
				if (utilityBinArraySU[i] != null
						&& utilityBinArraySU[i].sumsu >= lminUtil) {
					SU SUtemp = new SU();
					int maxIndex = -1;
					if (utilityBinArrayLU[newNamesToOldNames[i]].luPeriod
							.size() > 0) {
						for (Period luperiod : utilityBinArrayLU[newNamesToOldNames[i]].luPeriod) {
							int indexstart = luperiod.indexStart;
							int indexend = luperiod.indexEnd;
							if (maxIndex >= indexstart) {
								for (int p = maxIndex + 1; p <= indexend; p++) {
									SUtemp.pairSU
											.add(utilityBinArraySU[i].pairSU
													.get(p));
								}
								maxIndex = indexend;
							} else {
								for (int p = indexstart; p <= indexend; p++) {
									SUtemp.pairSU
											.add(utilityBinArraySU[i].pairSU
													.get(p));
								}
								maxIndex = indexend;
							}

						}
					}
					if (SUtemp.pairSU.size() > 0) {
						generatePeriodOfSU(SUtemp, lminUtil, minLen);
						if (SUtemp.suPeriod.size() > 0) {
							itemsToExplore.add(i);
						}
						if (SUtemp.uPeriod.size() > 0) {
							output(0, newNamesToOldNames[i],
									utilityBinArraySU[i].sumUtility,
									SUtemp.uPeriod);
						}
					}

				}

			}
		}

		// If in debug mode, show the list of promising items
		if (DEBUG) {
			System.out
					.println("------------------ itemsToExplore------------------- ");
			for (Integer integer : itemsToExplore) {
				System.out.println(newNamesToOldNames[integer]);
			}
			System.out.println(itemsToExplore.toString());
			System.out
					.println("---------------ts su suPeriod of single item-------------------");
			for (int i = 1; i < utilityBinArraySU.length; i++) {
				System.out.println("item : " + newNamesToOldNames[i] + " u : "
						+ utilityBinArraySU[i].sumUtility);
				System.out.print("suPeriod : ");
				for (Period period : utilityBinArraySU[i].suPeriod) {
					System.out.print("[" + period.start + "," + period.end
							+ "]");
				}
				System.out.println();
				System.out.print("uPeriod : ");
				for (Period period : utilityBinArraySU[i].uPeriod) {
					System.out.print("[" + period.start + "," + period.end
							+ "]");
				}
				System.out.println();
				for (PairSU pair : utilityBinArraySU[i].pairSU) {
					System.out.println("ts : " + pair.timestamp + " su : "
							+ pair.suOftimestamp + " u : " + pair.uOftimestamp);
				}
			}
		}

		// //======
		// Recursive call to the algorithm
		// If subtree utility pruning is activated
		if (activateSubtreeUtilityPruning) {

			// We call the recursive algorithm with the database, secondary
			// items and primary items
			backtrackingEFIM(dataset.getTransactions(), itemsToKeep,
					itemsToExplore, 1);
		} else {
			// We call the recursive algorithm with the database and secondary
			// items
			backtrackingEFIM(dataset.getTransactions(), itemsToKeep,
					itemsToKeep, 1);
		}

		// record the end time
		endTimestamp = System.currentTimeMillis();

		// close the output file
		if (writer != null) {
			writer.close();
		}

		// check the maximum memory usage
		MemoryLogger.getInstance().checkMemory();

		// return the set of high-utility itemsets
		return highUtilityItemsets;
	}

	/**
	 * Implementation of Insertion sort for sorting a list of items by
	 * increasing order of TWU. This has an average performance of O(n log n)
	 * 
	 * @param items
	 *            list of integers to be sorted
	 * @param items
	 *            list the utility-bin array indicating the TWU of each item.
	 */
	public static void insertionSort(List<Integer> items, LU[] itemToTWU) {
		// the following lines are simply a modified an insertion sort

		for (int j = 1; j < items.size(); j++) {
			Integer itemJ = items.get(j);
			int i = j - 1;
			Integer itemI = items.get(i);

			// we compare the twu of items i and j
			Integer comparison = itemToTWU[itemI].sumlu
					- itemToTWU[itemJ].sumlu;
			// if the twu is equal, we use the lexicographical order to decide
			// whether i is greater
			// than j or not.
			if (comparison == 0) {
				comparison = (itemI - itemJ);
			}

			while (comparison > 0) {
				items.set(i + 1, itemI);

				i--;
				if (i < 0) {
					break;
				}

				itemI = items.get(i);
				comparison = itemToTWU[itemI].sumlu - itemToTWU[itemJ].sumlu;
				// if the twu is equal, we use the lexicographical order to
				// decide whether i is greater
				// than j or not.
				if (comparison == 0) {
					comparison = itemI - itemJ;
				}
			}
			items.set(i + 1, itemJ);
		}
	}

	public void generatePeriodOfLU(LU luSet, int lminUtil, int minLen) {
		long initTime = System.currentTimeMillis();
		int lu = 0;
		int winEnd = 0;

		// these flags indicates if the first window is a PHUI period or LHUI
		// period
		boolean luPreflag = false;
		// find first window
		for (; winEnd < luSet.pairLU.size()
				&& luSet.pairLU.get(winEnd).timestamp < luSet.pairLU.get(0).timestamp
						+ minLen; winEnd++) {
			lu += luSet.pairLU.get(winEnd).luOftimestamp;
		}
		if (lu > lminUtil)
			luPreflag = true;

		// slide window to find all period information
		slideWindowOfLU(luSet, winEnd, lminUtil, lu, luPreflag, minLen);
		windowConsumerLU += (System.currentTimeMillis() - initTime);
	}

	private void slideWindowOfLU(LU luSet, int winEnd, int lminUtil, int lu,
			boolean luPreflag, int minLen) {
		int beginIndex = 0, endIndex = winEnd;
		for (int i = 0; i < luSet.pairLU.size();) {
			int x, y;

			for (y = i; y < luSet.pairLU.size()
					&& luSet.pairLU.get(y).timestamp == luSet.pairLU.get(i).timestamp; y++) {
				lu -= luSet.pairLU.get(y).luOftimestamp;
			}
			i = y;
			for (x = winEnd; x < luSet.pairLU.size()
					&& luSet.pairLU.get(x).timestamp < luSet.pairLU.get(y).timestamp
							+ minLen; x++) {
				lu += luSet.pairLU.get(x).luOftimestamp;
				winEnd = x + 1;
			}

			// add the high utility period that iUtil>minutil
			if (luPreflag) {
				if (lu < lminUtil) {

					luSet.luPeriod.add(new Period(
							luSet.pairLU.get(beginIndex).timestamp,
							luSet.pairLU.get(endIndex - 1).timestamp,
							beginIndex, endIndex - 1));
					luPreflag = false;
				} else
					endIndex = winEnd;
			} else {
				if (lu >= lminUtil) {
					luPreflag = true;
					beginIndex = i;
					endIndex = winEnd;
				}
			}
		}

	}

	public void generatePeriodOfSU(SU suSet, int lminUtil, int minLen) {
		long initTime = System.currentTimeMillis();
		int su = 0, u = 0;
		int winEnd = 0;

		// these flags indicates if the first window is a PHUI period or LHUI
		// period
		boolean suPreflag = false;
		boolean uPreflag = false;

		// find first window
		for (; winEnd < suSet.pairSU.size()
				&& suSet.pairSU.get(winEnd).timestamp < suSet.pairSU.get(0).timestamp
						+ minLen; winEnd++) {
			su += suSet.pairSU.get(winEnd).suOftimestamp;
			u += suSet.pairSU.get(winEnd).uOftimestamp;
		}

		if (su > lminUtil)
			suPreflag = true;
		if (u > lminUtil)
			uPreflag = true;
		// slide window to find all period information
		slideWindowOfSU(suSet, winEnd, lminUtil, su, u, suPreflag, uPreflag,
				minLen);
		windowConsumer += (System.currentTimeMillis() - initTime);
	}

	private void slideWindowOfSU(SU suSet, int winEnd, int lminUtil, int su,
			int u, boolean suPreflag, boolean uPreflag, int minLen) {
		int beginIndex = 0, endIndex = winEnd, uBeginIndex = 0, uEndIndex = winEnd;
		for (int i = 0; i < suSet.pairSU.size();) {
			int x, y;

			for (y = i; y < suSet.pairSU.size()
					&& suSet.pairSU.get(y).timestamp == suSet.pairSU.get(i).timestamp; y++) {
				su -= suSet.pairSU.get(y).suOftimestamp;
				u -= suSet.pairSU.get(y).uOftimestamp;
			}
			i = y;
			for (x = winEnd; x < suSet.pairSU.size()
					&& suSet.pairSU.get(x).timestamp < suSet.pairSU.get(y).timestamp
							+ minLen; x++) {
				su += suSet.pairSU.get(x).suOftimestamp;
				u += suSet.pairSU.get(x).uOftimestamp;
				winEnd = x + 1;
			}

			// add the high utility period that iUtil>minutil
			if (suPreflag) {
				if (su < lminUtil) {
					suSet.suPeriod.add(new Period(
							suSet.pairSU.get(beginIndex).timestamp,
							suSet.pairSU.get(endIndex - 1).timestamp));
					suPreflag = false;
				} else
					endIndex = winEnd;
			} else {
				if (su > lminUtil) {
					suPreflag = true;
					beginIndex = i;
					endIndex = winEnd;
				}
			}

			if (uPreflag) {
				if (u < lminUtil) {
					suSet.uPeriod.add(new Period(
							suSet.pairSU.get(uBeginIndex).timestamp,
							suSet.pairSU.get(uEndIndex - 1).timestamp));
					uPreflag = false;
				} else
					uEndIndex = winEnd;
			} else {
				if (u > lminUtil) {
					uPreflag = true;
					uBeginIndex = i;
					uEndIndex = winEnd;
				}
			}
		}

	}
	/**
	 * Recursive method to find all high-utility itemsets
	 * 
	 * @param the
	 *            list of transactions containing the current prefix P
	 * @param itemsToKeep
	 *            the list of secondary items in the p-projected database
	 * @param itemsToExplore
	 *            the list of primary items in the p-projected database
	 * @param the
	 *            current prefixLength
	 * @throws IOException
	 *             if error writing to output file
	 */
	private void backtrackingEFIM(List<Transaction> transactionsOfP,
			List<Integer> itemsToKeep, List<Integer> itemsToExplore,
			int prefixLength) throws IOException {
		
		// update the number of candidates explored so far
		candidateCount += itemsToExplore.size();
		// ======== for each frequent item e =============
		for (int j = 0; j < itemsToExplore.size(); j++) {
			Integer e = itemsToExplore.get(j);

			// ========== PERFORM INTERSECTION =====================
			// Calculate transactions containing P U {e}
			// At the same time project transactions to keep what appears after
			// "e"
			List<Transaction> transactionsPe = new ArrayList<Transaction>();

			// this variable is to record the time for performing intersection
			long timeFirstIntersection = System.currentTimeMillis();

			// For each transaction
			for (Transaction transaction : transactionsOfP) {
				
				// Increase the number of transa ction read
				transactionReadingCount++;

				// To record the time for performing binary searh
				long timeBinaryLocal = System.currentTimeMillis();

				// we remember the position where e appears.
				// we will call this position an "offset"
				int positionE = -1;
				// Variables low and high for binary search
				int low = transaction.offset;
				int high = transaction.items.length - 1;

				// perform binary search to find e in the transaction
				while (high >= low) {
					int middle = (low + high) >>> 1; // divide by 2
					if (transaction.items[middle] < e) {
						low = middle + 1;
					} else if (transaction.items[middle] == e) {
						positionE = middle;
						break;
					} else {
						high = middle - 1;
					}
				}
				// record the time spent for performing the binary search
				timeBinarySearch += System.currentTimeMillis()
						- timeBinaryLocal;

				// if 'e' was found in the transaction
				if (positionE > -1) {
					// Otherwise, if merging has been deactivated
					// then we just create the projected transaction
					Transaction projectedTransaction = new Transaction(
							transaction, positionE);
					// we put the projected transaction in the projected
					// database of Pe	
					transactionsPe.add(projectedTransaction);
					transaction.offset = positionE;
				} else {
					transaction.offset = low;
				}
			}
			
			// remember the total time for peforming the database projection
			timeIntersections += (System.currentTimeMillis() - timeFirstIntersection);

			// Append item "e" to P to obtain P U {e}
			// but at the same time translate from new name of "e" to its old
			// name
			temp[prefixLength] = newNamesToOldNames[e];

			// ==== Next, we will calculate the Local Utility and Sub-tree
			// utility of
			// all items that could be appended to PU{e} ====
			useUtilityBinArraysToCalculateUpperBounds(transactionsPe, j,
					itemsToKeep);
			// we now record time for identifying promising items
			long initialTime = System.currentTimeMillis();

			// We will create the new list of secondary items
			List<Integer> newItemsToKeep = new ArrayList<Integer>();
			// We will create the new list of primary items
			List<Integer> newItemsToExplore = new ArrayList<Integer>();
			
			if(!OPT)
			for (int k = j + 1; k < itemsToKeep.size(); k++) {
				int itemk = itemsToKeep.get(k);
				boolean signToSu = true;
				if (utilityBinArrayLU[itemk] != null
						&& utilityBinArrayLU[itemk].sumlu >= lminUtil) {
					generatePeriodOfLU(utilityBinArrayLU[itemk], lminUtil,
							minLen);
					if (utilityBinArrayLU[itemk].luPeriod.size() > 0) {
						newItemsToKeep.add(itemk);
					} else {
						signToSu = false;
					}

					if (signToSu && utilityBinArraySU[itemk] != null
							&& utilityBinArraySU[itemk].sumsu >= lminUtil) {
						SU SUtemp = new SU();
						int maxIndex = -1;
						if (utilityBinArrayLU[itemk].luPeriod.size() > 0) {
							for (Period luperiod : utilityBinArrayLU[itemk].luPeriod) {
								int indexstart = luperiod.indexStart;

								int indexend = luperiod.indexEnd;

								if (maxIndex >= indexstart) {
									for (int p = maxIndex + 1; p <= indexend; p++) {
										SUtemp.pairSU
												.add(utilityBinArraySU[itemk].pairSU
														.get(p));
									}
									maxIndex = indexend;
								} else {
									for (int p = indexstart; p <= indexend; p++) {
										SUtemp.pairSU
												.add(utilityBinArraySU[itemk].pairSU
														.get(p));
									}
									maxIndex = indexend;
								}

							}
						}

						if (SUtemp.pairSU.size() > 0) {
							generatePeriodOfSU(SUtemp, lminUtil, minLen);
							if (SUtemp.uPeriod.size() > 0)
								output(prefixLength, newNamesToOldNames[itemk],
										utilityBinArraySU[itemk].sumUtility,
										SUtemp.uPeriod);
							if (SUtemp.suPeriod.size() > 0) {
								newItemsToExplore.add(itemk);
							}
						}

					}
				}
			}
			else for (int k = j + 1; k < itemsToKeep.size(); k++) {
				int itemk = itemsToKeep.get(k);
				boolean signToLU = true;
				if (utilityBinArraySU[itemk] != null
						 && utilityBinArraySU[itemk].sumsu >= lminUtil) {
					generatePeriodOfSU(utilityBinArraySU[itemk], lminUtil,
							minLen);
					if (utilityBinArraySU[itemk].uPeriod.size() > 0)
						output(prefixLength, newNamesToOldNames[itemk],
								utilityBinArraySU[itemk].sumUtility,
								utilityBinArraySU[itemk].uPeriod);
					if (utilityBinArraySU[itemk].suPeriod.size() > 0) {
						signToLU = false;
						newItemsToExplore.add(itemk);
						newItemsToKeep.add(itemk);
					}
				}
				if (signToLU) {
					if (utilityBinArrayLU[itemk] != null
							&& utilityBinArrayLU[itemk].sumlu >= lminUtil) {
						generatePeriodOfLU(utilityBinArrayLU[itemk], lminUtil,
								minLen);
						if (utilityBinArrayLU[itemk].luPeriod.size() > 0)
							newItemsToKeep.add(itemk);
					}
				}
			}
			// update the total time for identifying promising items
			timeIdentifyPromisingItems += (System.currentTimeMillis() - initialTime);

			// === recursive call to explore larger itemsets
			if (activateSubtreeUtilityPruning) {
				// if sub-tree utility pruning is activated, we consider primary
				// and secondary items
				backtrackingEFIM(transactionsPe, newItemsToKeep,
						newItemsToExplore, prefixLength + 1);
			} else {
				// if sub-tree utility pruning is deactivated, we consider
				// secondary items also
				// as primary items
				backtrackingEFIM(transactionsPe, newItemsToKeep,
						newItemsToKeep, prefixLength + 1);
			}
		}

		// check the maximum memory usage for statistics purpose
		MemoryLogger.getInstance().checkMemory();
	}

	/**
	 * Scan the initial database to calculate the local utility of each item
	 * using a utility-bin array
	 * 
	 * @param dataset
	 *            the transaction database
	 */
	public void useUtilityBinArrayToCalculateLocalUtilityFirstTime(
			Dataset dataset) {

		// Initialize utility bins for all items
		utilityBinArrayLU = new LU[dataset.getMaxItem() + 1];

		// Scan the database to fill the utility bins
		// For each transaction
		for (Transaction transaction : dataset.getTransactions()) {
			int timestamp = transaction.timeStamp;
			// for each item
			for (Integer item : transaction.getItems()) {
				if (utilityBinArrayLU[item] == null) {
					utilityBinArrayLU[item] = new LU();
					PairLU pairLU = new PairLU(timestamp, 0);
					utilityBinArrayLU[item].pairLU.add(pairLU);
				}
				int previousIndex = utilityBinArrayLU[item].pairLU.size() - 1;
				if (utilityBinArrayLU[item].pairLU.get(previousIndex).timestamp == timestamp) {
					utilityBinArrayLU[item].pairLU.get(previousIndex).luOftimestamp = utilityBinArrayLU[item].pairLU
							.get(previousIndex).luOftimestamp
							+ transaction.transactionUtility;
					utilityBinArrayLU[item].sumlu += transaction.transactionUtility;
				} else {
					PairLU pairLU = new PairLU(timestamp,
							transaction.transactionUtility);
					utilityBinArrayLU[item].pairLU.add(pairLU);
					utilityBinArrayLU[item].sumlu += transaction.transactionUtility;
				}

			}
		}
	}

	/**
	 * Scan the initial database to calculate the sub-tree utility of each item
	 * using a utility-bin array
	 * 
	 * @param dataset
	 *            the transaction database
	 */
	public void useUtilityBinArrayToCalculateSubtreeUtilityFirstTime(
			Dataset dataset) {

		int sumSU;
		// Scan the database to fill the utility-bins of each item
		// For each transaction
		for (Transaction transaction : dataset.getTransactions()) {
			// We will scan the transaction backward. Thus,
			// the current sub-tree utility in that transaction is zero
			// for the last item of the transaction.
			sumSU = 0;
			int timestamp = transaction.timeStamp;
			// For each item when reading the transaction backward
			for (int i = transaction.getItems().length - 1; i >= 0; i--) {
				// get the item
				Integer item = transaction.getItems()[i];
				sumSU += transaction.getUtilities()[i];
				if (utilityBinArraySU[item] == null) {
					utilityBinArraySU[item] = new SU();
					PairSU pairSU = new PairSU(timestamp, 0, 0);
					utilityBinArraySU[item].pairSU.add(pairSU);
				}
				int previousIndex = utilityBinArraySU[item].pairSU.size() - 1;
				if (utilityBinArraySU[item].pairSU.get(previousIndex).timestamp == timestamp) {
					utilityBinArraySU[item].pairSU.get(previousIndex).suOftimestamp += sumSU;
					utilityBinArraySU[item].pairSU.get(previousIndex).uOftimestamp += transaction
							.getUtilities()[i];
					utilityBinArraySU[item].sumUtility += transaction
							.getUtilities()[i];
					utilityBinArraySU[item].sumsu += sumSU;
				} else {
					PairSU pairSU = new PairSU(timestamp, sumSU,
							transaction.getUtilities()[i]);
					utilityBinArraySU[item].pairSU.add(pairSU);
					utilityBinArraySU[item].sumUtility += transaction
							.getUtilities()[i];
					utilityBinArraySU[item].sumsu += sumSU;
				}

				// we add the current sub-tree utility to the utility-bin of the
				// item
				// utilityBinArraySU[item] += sumSU;
			}
		}
	}
	
	/**
	 * Utilize the utility-bin arrays to calculate the sub-tree utility and
	 * local utility of all items that can extend itemset P U {e}
	 * 
	 * @param transactions
	 *            the projected database for P U {e}
	 * @param j
	 *            the position of j in the list of promising items
	 * @param itemsToKeep
	 *            the list of promising items
	 */
	private void useUtilityBinArraysToCalculateUpperBounds(
			List<Transaction> transactionsPe, int j, List<Integer> itemsToKeep) {

		// we will record the time used by this method for statistics purpose
		long initialTime = System.currentTimeMillis();

		// For each promising item > e according to the total order
		for (int i = j + 1; i < itemsToKeep.size(); i++) {
			Integer item = itemsToKeep.get(i);
			// We reset the utility bins of that item for computing the sub-tree
			// utility and
			utilityBinArraySU[item] = null;
			utilityBinArrayLU[item] = null;
		}

		int sumRemainingUtility;
		// for each transaction
		for (Transaction transaction : transactionsPe) {
			// count the number of transactions read
			transactionReadingCount++;

			// We reset the sum of reamining utility to 0;
			sumRemainingUtility = 0;
			// we set high to the last promising item for doing the binary
			// search
			int high = itemsToKeep.size() - 1;

			// for each item in the transaction that is greater than i when
			// reading the transaction backward
			// Note: >= is correct here. It should not be >.
			for (int i = transaction.getItems().length - 1; i >= transaction.offset; i--) {
				// get the item
				int item = transaction.getItems()[i];
				// We will check if this item is promising using a binary search
				// over promising items.
				// This variable will be used as a flag to indicate that we
				// found the item or not using the binary search
				boolean contains = false;
				// we set "low" for the binary search to the first promising
				// item position
				int low = 0;
				// do the binary search
				while (high >= low) {
					int middle = (low + high) >>> 1; // divide by 2
					int itemMiddle = itemsToKeep.get(middle);
					if (itemMiddle == item) {
						// if we found the item, then we stop
						contains = true;
						break;
					} else if (itemMiddle < item) {
						low = middle + 1;
					} else {
						high = middle - 1;
					}
				}
				// if the item is promising
				if (contains) {
					// We add the utility of this item to the sum of remaining
					// utility
					int timestamp = transaction.timeStamp;
					sumRemainingUtility += transaction.getUtilities()[i];

					if (utilityBinArrayLU[item] == null) {
						utilityBinArrayLU[item] = new LU();
						PairLU pairLU = new PairLU(timestamp, 0);
						utilityBinArrayLU[item].pairLU.add(pairLU);
					}
					int previousIndex_LU = utilityBinArrayLU[item].pairLU
							.size() - 1;
					if (utilityBinArrayLU[item].pairLU.get(previousIndex_LU).timestamp == timestamp) {
						utilityBinArrayLU[item].pairLU.get(previousIndex_LU).luOftimestamp += transaction.transactionUtility
								+ transaction.prefixUtility;
						utilityBinArrayLU[item].sumlu += transaction.transactionUtility
								+ transaction.prefixUtility;
					} else {
						PairLU pairLU = new PairLU(timestamp,
								transaction.transactionUtility
										+ transaction.prefixUtility);
						utilityBinArrayLU[item].pairLU.add(pairLU);
						utilityBinArrayLU[item].sumlu += transaction.transactionUtility
								+ transaction.prefixUtility;
					}
					if (utilityBinArraySU[item] == null) {
						utilityBinArraySU[item] = new SU();
						PairSU pairSU = new PairSU(timestamp, 0, 0);
						utilityBinArraySU[item].pairSU.add(pairSU);
					}
					int previousIndex_SU = utilityBinArraySU[item].pairSU
							.size() - 1;
					if (utilityBinArraySU[item].pairSU.get(previousIndex_SU).timestamp == timestamp) {
						
						utilityBinArraySU[item].pairSU.get(previousIndex_SU).suOftimestamp += sumRemainingUtility
								+ transaction.prefixUtility;
						utilityBinArraySU[item].pairSU.get(previousIndex_SU).uOftimestamp += transaction
								.getUtilities()[i] + transaction.prefixUtility;
						utilityBinArraySU[item].sumUtility += transaction
								.getUtilities()[i] + transaction.prefixUtility;
						utilityBinArraySU[item].sumsu += sumRemainingUtility
								+ transaction.prefixUtility;
					} else {
						PairSU pairSU = new PairSU(
								timestamp,
								sumRemainingUtility + transaction.prefixUtility,
								transaction.getUtilities()[i]
										+ transaction.prefixUtility);
						utilityBinArraySU[item].pairSU.add(pairSU);
						utilityBinArraySU[item].sumUtility += transaction
								.getUtilities()[i] + transaction.prefixUtility;
						utilityBinArraySU[item].sumsu += sumRemainingUtility
								+ transaction.prefixUtility;
					}
					
				}
			}
		}
		// we update the time for database reduction for statistics purpose
		timeDatabaseReduction += (System.currentTimeMillis() - initialTime);
	}

	/**
	 * Save a high-utility itemset to file or memory depending on what the user
	 * chose.
	 * 
	 * @param itemset
	 *            the itemset
	 * @throws IOException
	 *             if error while writting to output file
	 */
	private void output(int tempPosition, int e, int utility,
			ArrayList<Period> period) throws IOException {
		patternCount++;

		// if user wants to save the results to memory
		if (writer == null) {
			// we copy the temporary buffer into a new int array
			int[] copy = new int[tempPosition + 1];
			System.arraycopy(temp, 0, copy, 0, tempPosition + 1);
			highUtilityItemsets.addItemset(new Itemset(copy, utility),
					copy.length);
		} else {
			StringBuffer buffer = new StringBuffer();
			for (int i = 1; i <= tempPosition; i++) {
				buffer.append(temp[i]);
				buffer.append(' ');
			}
			buffer.append(e);
			// append the utility of the itemset
			buffer.append(" #UTIL: ");
			buffer.append(utility);
			for (Period period2 : period) {
				buffer.append(" [" + period2.start + "," + period2.end + "]");
			}

			// write the stringbuffer to file and create a new line
			// so that we are ready for writing the next itemset.
			writer.write(buffer.toString());
			writer.newLine();
		}
	}

	/**
	 * Print statistics about the latest execution of the EFIM algorithm.
	 */
	public void printStats() {

		System.out.println("========== Fast-LHUI ============");
		System.out.println(" lminUtil ~ " + lminUtil);
		System.out.println(" minLen ~ " + minLen);
		System.out.println(" Local HUIs count ~ " + patternCount);
		System.out.println(" Total time ~ " + (endTimestamp - startTimestamp)
				/ 1000.0 + " s");
		// if in debug mode, we show more information
		if (DEBUG) {
			System.out.println(" Transaction merge count ~ " + mergeCount);
			System.out.println(" Transaction read count ~ "
					+ transactionReadingCount);

			System.out.println(" Time intersections ~ " + timeIntersections
					+ " ms");
			System.out.println(" Time database reduction ~ "
					+ timeDatabaseReduction + " ms");
			System.out.println(" Time promising items ~ "
					+ timeIdentifyPromisingItems + " ms");
			System.out.println(" Time binary search ~ " + timeBinarySearch
					+ " ms");
			System.out.println(" Time sort ~ " + timeSort + " ms");
		}
		System.out.println(" Max memory:"
				+ MemoryLogger.getInstance().getMaxMemory());
		System.out.println(" Candidate count : " + candidateCount);
		System.out.println(" Time binary search ~ " + timeBinarySearch/1000.0
				+ " s");
		System.out.println(" timeDatabaseReduction : " + timeDatabaseReduction/1000.0+ " s");
		System.out.println("=====================================");
	}
}
