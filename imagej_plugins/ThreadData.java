

public class ThreadData
{
	public static int use_threads = 0;
	String command = null;
	int relaunches = 0;
	int job = -1;
	int rank = -1;
	int worker = -1; // counting the entries in nodes.txt
	long zipImgLen = -1;
	double performance = -1.0;
	double load = -1.0;
	boolean can_write = false;
}
