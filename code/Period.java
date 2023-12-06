package algo;


public class Period {
	public int start;
	public int end;
	public Integer indexStart;
	public Integer indexEnd;
	public Period(int start, int end) {
		super();
		this.start = start;
		this.end = end;
	}
	@Override
	public String toString() {
		return "[" + start + ", " + end + "]";
	}
	public Period(int start, int end, Integer indexStart, Integer indexEnd) {
		super();
		this.start = start;
		this.end = end;
		this.indexStart = indexStart;
		this.indexEnd = indexEnd;
	}
}

