import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.FileSaver;
import ij.text.TextWindow;
import java.awt.*;
import java.awt.geom.*;
import ij.plugin.*;
import java.util.Random;
import java.util.Scanner;
import java.util.Arrays;
import java.util.Vector;
import java.util.Locale;
import java.util.HashMap;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.IOException;

import java.nio.ByteBuffer;
import java.text.DecimalFormat;

import java.util.Date;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

import krami.SharedMerit;


public class Hex_Magic implements PlugIn {

    //Global variables (as few as possible)
    long seed = System.currentTimeMillis();
    Random random = null; // not to be used until seed is confirmed by user
	//Permutator permutator = null; //not be used before modelsize is known
	Data_Generator datagen = null; //in case we want to make our own data
	static byte[] bytes = new byte[8]; //for FAST stream I/O operations (at least one double must fit)
	ByteBuffer buffer = ByteBuffer.wrap(bytes);
    static int mutex = 0;
    int[][] permutation = null;
    
    public static java.text.DecimalFormat format = new DecimalFormat("0.000000");
    DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	Date date = new Date();
	String sourcepath = IJ.getDirectory("imagej"); //default is home of ImageJ
    String[] Modii = {"run", "write job files", "read job files","models"};
    String Modus = null; //one of the modii
    
    public static SharedMerit sharedmerit = null;
   
    
    double brightness_offset = 3; //barely sufficient for void in 8 ring
    double datamean0 = 1;
    double data_bg0 = 0.0;
    double data_bg_std0 = 0.0;
    long datasum0;
    int dataStSize0 = 0;
    long[] bigSumImg = null;
	double modelmin0;
	double modelmax0;
	double modelmean0;
	double noise_level = 0.0;
	double clone_weight = 0.05;
    
    double default_weight = 0.015;
    double init_lattice_weight = 0.75;
    double[] modelweight, modelweightupd, sym_factors;
    //double sym_fac_exponent = 1.0;
    int[] cases;
    //final String sparseclient = "_sparse";
    
    String[] modelStatus = {"Continue","Finished","Fixed","Lattice","Reset", "Noisify", "Error"}; //"" is shorthand for "Continue"
    enum ModelState {CONTINUE, FINISHED, FIXED, LATTICE, RESET, NOISIFY, NOMODELSTATE;};
    ModelState[] dead; //index of modelState
    //the state of a model decides what can happen
    //"Finished" only weight adjustments, no more pixel tests
    //"Reset" it to copy of Lattice, add noise and "Continue"
    //"Noisify" introduce some noise and then "Continue"
    //"Continue" keep ajusting pixels and weights, "Finished" only adjust weights
    //"Lattice" no further weight or pixel changes
    int defaultmodelnum = 4;
    int modelnum = defaultmodelnum ;  //number of model images
    int histonum = 0;
    int modelGraylevels = 128;
    double modelBrightness = -4.0;

    enum Strategy {PIXELADJ, WEIGHTADJ, PIXEL_EM, MODEL_EM, TOTAL_EM, PIX_MOVES, NOSTRATEGY, CONFIRMED_PIXEL,POLL_HISTO};
    Strategy strategy = Strategy.PIXELADJ;
    //Binary optimization of Single Pixel
	int initthreadnum = 0;

	int renew = 0; //signal to pick new strategy
	int cmsok = 0; //a strategy was picked


    int wrapnum = 0;
    int pause_at_wrap = -1;
    int new_wraps = 0;
    int init_target_wraps = -1;
    int target_wraps = init_target_wraps; //Quit after so many wraps have been completed
    int init_target_models = -1;
    int target_models = init_target_models; //Quit after so many or models have finished
   
	int num_pixels_active = 0;
	int num_pixels_total = 0;
	int weight_trials = 0;
	int move_trials = 0;
	int the_slice = 1;
	int multiActivations = 1;
	int multiActivations0 = 1;
	int wobble = 0;

    int[] unsc;
    int unscwr = 0;

	double initpvaloffset = 0.0;
	double pvaloffset = initpvaloffset;
	double pvalscaling = 1.0;
	double noisification_level = 0.1;
	boolean spiral = false;
	boolean append_Pval = false;

	int spiralDelta = 1; //1 .. growing, .. -1 shrinking
	int spiraldirection = 0; //0 .. left, 1 .. down, 2 .. right, 3 .. up
	int spiralcount = -1; //Size of the spiral
	int spiralcountmax = 1; //edge of length, increases after direction 3 has reached max

	public static int bondlength = 4; //pixels per hexagon
	int initblocksize = 1;
	int blocksize = initblocksize;
	int blockarea = blocksize * blocksize;
	int new_blocksize = blocksize;
	int initblockdepth = 0;
	int blockdepth = initblockdepth;
	int new_blockdepth = blockdepth;
	int blockvolume = blockarea*blockdepth;
	int total_subframes = 0;
	int shadownum = 1;
							

	int blockid = -1; //used for accessing shuffled blocks
	//internal coordinates of the block
	int blockx = 0;
	int blocky = 0;
	int blockm = 0;
	//home coordinates of the datablock
	int homex = 0;
	int homey = 0;
	int homem = 0;

	public static int modelsize = -1;  // size of model image (width and height)
	int modelactiveinit = -120; //half the diameter corresponds to biggest possible distance
	int modelactive = modelactiveinit;
	int beamsize = -1;
	int beamarea = -1;


	int ongoing_Adj = -1;
	int incrs[]; // 0 .. fresh test, 1 .. increasing steps, 2 .. decreasing step
       
	int testval[];    //the current guess
	int oldvalinit[]; //the starting point
	
	int oldval[]; //the last guess
    int bestval[];//the best guess so far
    int secondbestval[]; //the secondbest guess so far
    int thirdbestval[]; //the third best guess so far
    int avg_blockval[];
    double Pval = 0.0;
    double highestPval = 0.0;
    double highestPval_at_blockinit = Double.NaN;

    double bestPval[]; //the highest Pval
    double secondbestPval[]; //the secondhighest Pval
	double thirdbestPval[]; //the third highest Pval

	double oldwght;
	double newwght;
	double bestwght;
	double secondbestwght;
	double thirdbestwght;
	double bestPvalW;
	double secondbestPvalW;
	double thirdbestPvalW;
	double wstep;
	int weightupdnum = 0;
	int weight_incrs;

	boolean external_merit = false;
	boolean fresh_symmetries = false;
    boolean writehist = false;
    boolean request_cache_report = false;
    boolean request_autosave = false;
    boolean pending_keep = false;
    boolean batch_run = false;
    boolean new_models = false;
    boolean noisy_init_models = false;
    boolean another_job = false;
    boolean exclusive_autosave = true;
    boolean keep_Idata = true;
    boolean real_log2 = false;
    boolean check_free_nodes = true;
    boolean clone_models = false;
    boolean skip_first_clone = true;
    boolean interactive_cloning_session = true;
    boolean blobs = false;
    boolean sym_scaling = false;
    boolean median_filtering = false;
    boolean do_full_EM = false;
    boolean do_pix_EM = false;
    boolean force_pix_EM = false;
    boolean collect_Histos = false;
    boolean use_beamprofile = false;
    boolean testing_avg = false;
    boolean reset_activity = false;
    boolean aggressive_weights = true;
    boolean lattice_hopping = true;
    //boolean sparse_data = false;
    public static int job_id = 0;
    public static int rank_id = 0;
    public static boolean duplicate_Idata = true;
    public static boolean skip_stats = false;
    int default_job = 0;

    int minfornextm = 0;
    int deadModelCount = 0;

	int dbg_weights = 0;
    int dbg_bisect = 0;
    int dbg_quadfit = 0;
    int dbg_moves = 0;
    int dbg_good_moves = 0;
    int pixel_tests = 0;

    int threadadd = 0;
    int max_threads = 0;
    int use_threads = 0;
    int skip_threads = 0;
    int streamnum = 0;
    int looptotal = 0;
    //int smooth_passes = 1;

    long timeoffset = 0;
    //double weightfr = 0.025;
    double rotfr = 0.0;
    double fitfr = 0.0;


	boolean autoptableoffset = true;
	boolean initbenchmark = false;
	boolean shuffle_frames = true;
	boolean shuffle_update = false;
	boolean blockdepth_priority = true;
	boolean flip_block = false;
	boolean alternating_scan_direction = false;
	boolean flip_model = false;
    boolean reverse_block_depth = false;
    boolean correlated_optimization = false;
    boolean quick_optimization = false;
    boolean rough_optimization = false;
    boolean enable_quad_fits  = true;
    boolean noisify_at_restart = false;
    boolean fixed_lattice_weight = false;
    boolean frozen_lattice = false;
    boolean pixel_done = false;
    boolean normalize_model_intensities = true;
    boolean previous_seed = false;
    boolean display_options = false;
    boolean fair_weights = false;
    boolean use_weights = true;
    boolean const_weights = false;
    boolean next_use_weights = true;
    boolean tell_use_weights = false;
    boolean init_weights = false;
    boolean use_max = false;
    boolean next_use_max = false;
    boolean tell_use_max = false;
    boolean black_init = false;
    boolean can_do_stats = false; //if Significance.txt should be written
	boolean master_shows_dialog = true;
	boolean suppress_further_dialog = false;
	
	public static boolean[] hex_mask = null;

	//boolean force_wrap = false;
    int debug_level = 2;
/*
	0 ... speed, 1 ... write Pval.txt, 2 ... autosave images every wrap, 3 ... also request Cachereport every wrap, 4 ... also write Job.bin,
	conversions from old options.txt
	false -> 0, true -> 1
*/
    

    int newdiffval=0;
    int reply = -1;
    int init_use_frames = 0;
    int use_frames = init_use_frames;

    double area_noise_level = 0.1;
	String oldbeamProfileTitle = "<none>";
	String beamProfileTitle; 
    String oldModeltitle = "";
    String Modeltitle;
    String oldDatatitle = "";
    String Datatitle;
    String oldPtableTitle = "";
    String PtableTitle;
    String master_PtableTitle; 
    String Weighttitle = "Weights.txt";
	String s = ""; //the current reply of a worker
	String scan_pattern = "noise";
	String user = System.getProperty("user.name");
	String prefix = "";
	String master_prefix = "";
	String outpath = "";
	
    int report_interval = 10000;
    int report_num = 3;
    int report_cols = 1;
    int weight_cols = 4;

	int thread_offset = 0;
	int MAX_THREADS = 255;
    
    //coords of active Hex
	//They are static so that other classes can use the test for coordinates
	public static int cq = 0;
	public static int cr = 0;
	public static int cx = cq - cr / 2;
	public static int cz = cr;
	public static int cy = -cx - cz;

    Vector<JobData> jobVec = new Vector<JobData>();   //job related data
    Vector<ThreadData> threadVec = new Vector<ThreadData>(); //thread related data
    Vector<String> batchpaths = new Vector<String>();
	Vector<String> prefix_per_path =  new Vector<String>();
	
	
	ImagePlus Idata = null; // original ImagePlus Stack
	ImageStack IdataSt = null;

	ImagePlus Diffimg = null;
	ImageStack DiffSt = null;

	ImagePlus nextModel = null;
	ImagePlus Models = null; // ImagePlus for model snapshots
	ImageStack ModelsSt = null;
	
	ImagePlus ext_Model = null;
	ImageStack ext_ModelSt = null;

	ImagePlus Ptable = null;
	ImageProcessor Ptablep = null;
	
	ImagePlus NewPtable = null;
	ImageProcessor NewPtablep = null;
	
	ImagePlus Updates = null;
	ImageStack UpdatesSt = null;
	short[] update_pixels = null;
	
	ImagePlus Activity = null;
	ImageStack ActivitySt = null;
	short[] activity_pixels = null;
	
	ImagePlus beamProfile = null;
	ImagePlus Histimg = null;
	
	
	
	//Diff_Matcher matcher = null;
	Diff_Matcher nextMatcher = null;
	
    private void resetGlobals()
    { //reset the needed globals for the batchrun
		sourcepath = IJ.getDirectory("imagej");
		Weighttitle = "Weights.txt";
		PtableTitle = null;
		wrapnum = 0;
		pause_at_wrap = -1;
		new_wraps = 0;
		Pval = 0.0;
		highestPval = 0.0;
		writehist = false;
		request_cache_report = false;
		pending_keep = false;
		new_models = false;
		normalize_model_intensities = false;
		master_shows_dialog = false; // has to be explicitely set in options.txt
		sym_scaling = false;
		median_filtering = false;
		do_full_EM = false;
		do_pix_EM = false;
		force_pix_EM = false;
		external_merit = false;
		
		collect_Histos = false;
		use_beamprofile = false;
		histonum = 0;
		total_subframes = 0;
		deadModelCount = 0;
		use_threads = 0;
		jobVec.clear();
		job_id = 0;
		rank_id = 0;
		dbg_weights = 0;
		dbg_bisect = 0;
		dbg_quadfit = 0;
		dbg_moves = 0;
		dbg_good_moves = 0;
		pixel_tests = 0;
		num_pixels_active = 0;
		num_pixels_total = 0;
		weight_trials = 0;
		move_trials = 0;
		looptotal = 0;
		timeoffset = 0;
		//internal coordinates of the block
		blockx = 0;
		blocky = 0;
		blockm = 0;
		//home coordinates of the datablock
		homex = 0;
		homey = 0;
		homem = 0;
		
	}
    
    public static void setHexcenter(int q, int r)
	{
		//in principle negative centers are allowed
		cq = q;
		cr = r;
		cx = q - r / 2;
		cz = r;
		cy = -cx - cz;
	}

	public static boolean hexagonaltrap(int[] cube, int N)
	{
		int lx = cube[0];
		int ly = cube[1];
		int lz = cube[2];
		boolean at_home = true;
		//apply periodic boundary conditions
		int runs = 0;
		do
		{
			at_home = true;
			if(lx >= N)
			{
				at_home = false;
				lx -= (2*N);
				ly -= (-N);
				lz -= (-N);
			}
			else if (lx < -N)
			{
				at_home = false;
				lx += (2*N);
				ly += (-N);
				lz += (-N);
			}

			if(ly > N)
			{
				at_home = false;
				lx -= (-N);
				ly -= (2*N);
				lz -= (-N);
			}
			else if (ly <= -N)
			{
				at_home = false;
				lx += (-N);
				ly += (2*N);
				lz += (-N);
			}

			if(lz >= N)
			{
				lx -= (-N);
				ly -= (-N);
				lz -= (2*N);
			}
			else if (lz < -N)
			{
				lx += (-N);
				ly += (-N);
				lz += (2*N);
			}
			++runs;
		} while (!at_home);

		cube[0] = lx;
		cube[1] = ly;
		cube[2] = lz;
		return (runs == 1);
	}

    public static boolean iswithinHex(int q, int r, int cradius)
	{
		if(q < 0 || r < 0)
		{	IJ.error("Hex_Magic.iswithinHex: image coordinates cannot be negative");}
		//cube coords of target hex
		int px = q - r / 2;
		int pz = r;
		int py = -px - pz;
		//offset hex in cube coords
		int dx = px - cx;
		int dy = py - cy;
		int dz = pz - cz;
		return (dx >= -cradius && dx <  cradius &&
				dy >  -cradius && dy <= cradius &&
				dz >= -cradius && dz <  cradius );
	}

	private static int distToCenter(int q, int r)
	{
		//cube coords of target hex
		int px = q - r / 2;
		int pz = r;
		int py = -px - pz;
		//offset hex in cube coords
		int dx = px - cx;
		int dy = py - cy;
		int dz = pz - cz;
		//dist to central position
		return ( ( Math.abs(dx) + Math.abs(dy) + Math.abs(dz) ) / 2);
	}


    public void run(String arg)
    {
        IJ.register(getClass());
        if(mutex == 1)
		{	
			IJ.error("Sorry, there is aleady another instance of" +  getClass().getSimpleName()  + " running!");
			return;
		}
        
        Locale.setDefault(Locale.UK);
        Datatitle = "";
		File sendjob = null;
		File nuke = null;
		File jobrun = null;
		File jobStats = null;
		File options = null;
		File batchfile = null;
		File nodes = null;
		File usenodes = null;
		Scanner usenodes_sc = null;

		java.io.FileOutputStream Pvalstream = null;
		java.io.PrintWriter Pvalhistory = null;

		java.io.FileOutputStream sendjobstream = null;
		java.io.PrintWriter sendjobwriter = null;
		
		java.io.FileOutputStream usenodes_stream = null;
		java.io.PrintWriter usenodes_writer = null;
		
		java.io.FileOutputStream nukestream = null;
		java.io.PrintWriter nukewriter = null;

        int impHeight;  //size of input image (pixel width or height)
        int impWidth;
        int histwidth = 1; 
		int histheight = 1;
        int stackSize;  //number of data images
        double modelmin = 0.0;
		double modelmax = -1.0;
		double datamean, dataStdDev, data_bg, data_bg_std;
		int  datamin, datamax;
		long datasum, data_num_hexes;
		String m_Datatitle = Datatitle;
		long m_seed = seed;
		String m_sourcepath = sourcepath;
		int m_init_use_frames = init_use_frames;
		int m_initthreadnum = initthreadnum;
		boolean m_previous_seed = previous_seed;
		TextWindow statusdisplay = null;
		long[] dataSumImg = null;
		//long[] bigSumImg  = null;
		boolean job_stats_ok = false;
		boolean pause_at_next_wrap = false;
		
		int threadnum;
		String[] pattern = {"spiral", "noise", "lines"};
		date.setTime(System.currentTimeMillis());
		IJ.log(getClass().getSimpleName() + " " + dateFormat.format(date) + "  " + user);
		
		String path = IJ.getDirectory(getClass().getSimpleName() + " Job0");
		String path0 = path;
		String master_path = path; 
		
		if(path == null) {	return;}
        
        
        IJ.log(path);
		// Get runtime
		java.lang.Runtime rt = java.lang.Runtime.getRuntime();
		//arrays for pid etc.
		java.lang.Process[] p = null;//new java.lang.Process[max_threads];  //pid
		//the process' output (our input stream- here we read responses from the worker thread)
		java.io.InputStream[] is = null;//new java.io.InputStream[max_threads];
		java.io.BufferedReader[] reader = null;//new java.io.BufferedReader[max_threads];
		//the process' input stream (our output - here we send data and commands to it)
		java.io.OutputStream[] os = null;//new java.io.OutputStream[max_threads + 1];
		java.io.PrintWriter[] pwos = null;//new java.io.PrintWriter[max_threads + 1];
		boolean[] can_write = null;//new boolean[max_threads + 1];
		batchfile = new File(path0 + "batchrun.txt");
		if(batchfile.isFile() && batchfile.setReadable(true))
		{
			boolean paths_ok =	parse_batchfile(0, batchfile);
			if(!paths_ok)
			{	return;}
			if(!batchpaths.isEmpty())
			{	
				batch_run = true;
				statusdisplay = new TextWindow(getClass().getSimpleName() + " status","will run in batch mode as long as " + batchfile.getName() + " is open",600,250);
			}
		}
		if(batchpaths.isEmpty())
		{
			IJ.log("empty batchpaths looking for options.txt");
			batchpaths.add(path0 + "options.txt");
			prefix_per_path.add(prefix);
			batch_run = false;
		}
		options = new File(batchpaths.get(0));
		prefix = prefix_per_path.get(0); //ok not really needed but somehow cleaner
		for(int batch_id = 0; batch_id < batchpaths.size(); ++batch_id)
		{		
			//if(batch_id > 0) //Will this prevent wrong path in a fresh options.txt
			resetGlobals();
			if( options.isFile() && options.setReadable(true) ) 
			{
				if(batch_id == 0) //do that early so that loadNodes knows about scanning the network
				{
					if(!loadOptions(options, true )) //the very first is always master
					{	return;}
					IJ.log("loaded " + options.getPath());
				}
			}
			else if(batch_run)
			{
				IJ.log("could not open: " + options.getPath());
				return;
			}
			
			boolean domainloop = true;
			boolean recent_save = false;
			boolean recent_img_save = false;
			
			max_threads = loadNodes(path0, true);
			if(max_threads < 1)
			{	
				IJ.log("There are no available workes. Please check for zombies");
				return;
			}
			if(max_threads > MAX_THREADS)
			{
				IJ.log("The limit of supported workers is " + MAX_THREADS + " . You may set it in options.txt");
			}
			p = new java.lang.Process[MAX_THREADS];  //pid
			//the process' output (our input stream- here we read responses from the worker thread)
			is = new java.io.InputStream[MAX_THREADS];
			reader = new java.io.BufferedReader[MAX_THREADS];
			//the process' input stream (our output - here we send data and commands to it)
			os = new java.io.OutputStream[MAX_THREADS + 1];
			pwos = new java.io.PrintWriter[MAX_THREADS + 1];
			can_write = new boolean[MAX_THREADS + 1];
						
			if(batch_run)
			{
				statusdisplay.append(batchpaths.get(batch_id) + " (" +  ( batch_id + 1 ) + " of " + batchpaths.size() + ")" );
				System.out.println("Running job " + ( batch_id + 1 ) + " of " + batchpaths.size() );
				IJ.log("Running job " + (batch_id + 1) + " of " + batchpaths.size());
			}
			try
			{
				do //another_job
				{
					boolean master = (job_id == 0);
					if(!master) {	another_job = false;}//dont add further jobs unless explicitely required 
					if(master || batch_run) 
					{	
						path = batchpaths.get(batch_id);
						prefix = prefix_per_path.get(batch_id);
					}
					else //neither master nor batchrun aka manual adding of another job
					{   
						path = IJ.getDirectory("Hex Magic Job" + job_id);
						if(path == null) break; // leave do another_job
					}
					
					jobVec.add(new JobData() );
					int firstthread = use_threads;
					if(!master) //append a new path rather than another options file
					{
						System.out.println("appending " + path );
						IJ.log("appending " + path );
						if(batch_run)
						{
							statusdisplay.append("appending " + path + " (" + ( batch_id + 1 ) + " of " + batchpaths.size() + ")");
						}
					}
					options = new File(path);
					String parentpath = options.getParent() + File.separator;
					jobVec.get(job_id).path = parentpath;
					jobrun = new File(parentpath + "jobrun.txt");
					sendjob = new File(parentpath + "sendjobs.sh");
					if(master)
					{
						if(outpath.equals(""))
						{
							master_path = parentpath;
						}
						else
						{
							master_path = outpath;
							outpath = ""; //outpath should not be sticky for sequential batchruns
						}
						nuke = new File(master_path + "killworkers.sh");
						jobStats = new File(master_path + "jobStats.txt");
						System.out.println("output = " + master_path);
					}
					
					
					
					if(options.isFile())
					{
						if(batch_id > 0)
						{
							if(!loadOptions(options, master ))
							{	return;}
							IJ.log("loaded " + options.getPath());
						}
					}
					else
					{
						IJ.log("missing file: " + options.getPath());
						//IJ.error("terminating batch mode, you may still manually launch the jobs");
						batch_run = false;
						sourcepath = options.getParentFile().getParent() + File.separator;
						IJ.log("assuming source path: " + sourcepath);
					}
					System.out.println("input  = " + sourcepath);
					if( (!( oldDatatitle.equals("<none>") || oldDatatitle.equals("<generate>") ) 
						&& !oldDatatitle.equals("") ) && (WindowManager.getImage(oldDatatitle) == null) )
					{
						File tmp = new File(sourcepath + oldDatatitle + ".zip");
						if(tmp.isFile())
						{
							Idata = new ImagePlus(tmp.getPath());
							IJ.log("opened " + tmp.getPath());
							Idata.show();
						}
						else
						{
							//IJ.log("could not open " + tmp.getPath());
							tmp = new File(sourcepath + oldDatatitle);
							if(tmp.isFile())
							{
								Idata = IJ.openVirtual(tmp.getPath());
								IJ.log("opened " + tmp.getPath());
								Idata.show();
							}
							else
							{
								IJ.log("could not open " + tmp.getPath());
							}
						}
					}
					if(WindowManager.getImage(oldModeltitle) == null)
					{
						File tmp = new File(parentpath + oldModeltitle);
						if(tmp.isFile())
						{
							Models = new ImagePlus(tmp.getPath());
							IJ.log("opened " + tmp.getPath());
							Models.show();
						}
						else
						{
							IJ.log("could not find " + tmp.getPath());
							tmp = new File(sourcepath + oldModeltitle);
							if(tmp.isFile())
							{
								Models = new ImagePlus(tmp.getPath());
								IJ.log("opened " + tmp.getPath());
								Models.show();
							}
						}
					}
					if((!oldbeamProfileTitle.equals("<none>")) && 
					    WindowManager.getImage(oldbeamProfileTitle) == null)
					{
						File tmp = new File(parentpath + oldbeamProfileTitle);
						if(tmp.isFile())
						{
							beamProfile = new ImagePlus(tmp.getPath());
							IJ.log("opened " + tmp.getPath());
							beamProfile.show();
						}
						else
						{
							IJ.log("could not find " + tmp.getPath());
							tmp = new File(sourcepath + oldbeamProfileTitle);
							if(tmp.isFile())
							{
								beamProfile = new ImagePlus(tmp.getPath());
								IJ.log("opened " + tmp.getPath());
								beamProfile.show();
							}
						}
					}
					if(WindowManager.getImage(oldPtableTitle) == null)
					{
						File tmp = new File(parentpath + oldPtableTitle);
						if(tmp.isFile())
						{
							Ptable = new ImagePlus(tmp.getPath());
							IJ.log("opened " + tmp.getPath());
							Thread.sleep(10);
							Ptable.show();
						}
						else
						{
							tmp = new File(sourcepath + oldPtableTitle);
							if(tmp.isFile())
							{
								Ptable = new ImagePlus(tmp.getPath());
								IJ.log("opened " + tmp.getPath());
								Thread.sleep(10);
								Ptable.show();
							}
						}
					}
			  
			  
					if(	(WindowManager.getImage(oldDatatitle) != null) 
						|| oldDatatitle.equals("<none>")
						|| oldDatatitle.equals("<generate>") )
					{	Datatitle = oldDatatitle;}
					else
					{
						if(batch_run)
						{
							IJ.error("failed to open Idata during batchrun\n" + oldDatatitle);
							String newdata = IJ.getFilePath("open Idata");
							if(newdata == null)
							{	return;}
							
							if(newdata.endsWith(".zip"))
							{
								Idata = new ImagePlus(newdata);
							}
							else
							{
								Idata = IJ.openVirtual(newdata);
							}
							
							
							
							Idata.show();
						}
					}
					if(WindowManager.getImage(oldModeltitle) != null)
					{	Modeltitle = oldModeltitle;	}
					else
					{	Modeltitle = "avg. data";}
					if(WindowManager.getImage(oldbeamProfileTitle) != null)
					{	beamProfileTitle = oldbeamProfileTitle;	}
					else
					{	beamProfileTitle = "<none>";}
					
					
					if(WindowManager.getImage(oldPtableTitle) != null)
					{	PtableTitle = oldPtableTitle;}
					else
					{	PtableTitle = "<new>";}
					
					IJ.log("Available workers: " + (threadadd - use_threads) );
					NonBlockingGenericDialog gd = new NonBlockingGenericDialog( getClass().getSimpleName() + " for " + user);
					gd.setSmartRecording(true);
					int[] idArray = WindowManager.getIDList();
					int idlen = 0;
					if(idArray != null)
					{	idlen = idArray.length;}
					String[] titleArray = new String[idlen + 2];
					for (int i = 0; i < idlen; ++i)
					{	titleArray[i] = WindowManager.getImage(idArray[i]).getTitle();}
					titleArray[idlen] = "<none>";
					titleArray[idlen + 1] = "<generate>";
					if(Datatitle.equals(""))
					{	Datatitle = titleArray[0];} //should be the oldest image
					
					String[] titleArrayPlus = new String[idlen + 2];
					titleArrayPlus[0] =	"avg. data";
					titleArrayPlus[1] = "none/noise";
					for (int i = 0; i < idlen; ++i)
					{	titleArrayPlus[i+2] = WindowManager.getImage(idArray[i]).getTitle();}
					
					String[] ptableArray = new String[idlen + 2];
					ptableArray[0] = "<new>";
					ptableArray[idlen+1] = "<none>";
					for (int i = 0; i < idlen; ++i)
					{	ptableArray[i+1] = WindowManager.getImage(idArray[i]).getTitle();}
					if(PtableTitle.equals(""))
					{	PtableTitle = ptableArray[0];	}
					
					String[] profileArray = new String[idlen + 1];
					profileArray[0] = "<none>";
					for (int i = 0; i < idlen; ++i)
					{	profileArray[i+1] = WindowManager.getImage(idArray[i]).getTitle();}
					if(WindowManager.getImage(beamProfileTitle) == null)
					{	beamProfileTitle = profileArray[0];}
					
					
					
					gd.addMessage(path);
					if(master) 
					{
						gd.addNumericField("Default number of models(<1 stacksize): ", modelnum, 0);
						gd.addNumericField("Rotated copies within pi/3 (>=1): ", shadownum, 0);
						gd.addNumericField("paralell model activations: ", multiActivations0, 0);
						gd.addNumericField("Model init. active diameter (n<0 -> -60*size/n): ", modelactiveinit, 0);
						gd.addNumericField("Model gray levels: ", modelGraylevels, 0);
						gd.addNumericField("Equivalent max data in model1(x<0 x*mean)", modelBrightness, 3);
						gd.addNumericField("Bond Length:", bondlength, 0);
						gd.addNumericField("Wobble radius: ",wobble,0);
						gd.addNumericField("Block size (<0 => x<->y):", initblocksize, 0);
						gd.addNumericField("block depth (<0 => reverse)", initblockdepth,0);
						gd.addNumericField("skip the first N available workers :", skip_threads, 0);
						gd.addNumericField("relative weight for clones :", clone_weight,3);
					}
					
					gd.addNumericField("Number of workers (<1: all):", initthreadnum, 0);
					if(master)
					{
						gd.addNumericField("worker ID offset :", thread_offset, 0);
						
						gd.addNumericField("Stop after nth stable modelds (-1,last)", init_target_models, 0);
						gd.addNumericField("Stop after nth wrap (-1,never)", init_target_wraps, 0);
					}
					gd.addNumericField("Number of sub frames to use (0,All)", init_use_frames, 0);
					if(master)
					{
						gd.addNumericField("initial model1 weight (w<0 => -w*fair)", init_lattice_weight, 3);
						gd.addCheckbox("Quit after first Pval", initbenchmark);
						if(batch_run)
						{	gd.addCheckbox("suppress all further dialogs", suppress_further_dialog);}
					}
					gd.addCheckbox("Shuffle frames", shuffle_frames);
					gd.addCheckbox("reuse previous random seed", previous_seed);
					
					String[] labels1 = {
						"alternating scan orientation",
						"rescale models",
						"blobs",
						"correlated optimization",
						
						"quick optimization",
						"fixed model1 weight",
						"frozen model1",
						"use weights",
						
						"fair weights",
						"const weights",
						"use max",
						"black initalization",
						
						"exclusive autosave",
						"keep Idata",
						"actual log2",
						"clone models"};
					boolean[] states1 = {   
						alternating_scan_direction,
						normalize_model_intensities,
						blobs,
						correlated_optimization,
						
						quick_optimization,
						fixed_lattice_weight, 
						frozen_lattice,
						use_weights,
						
						fair_weights,
						const_weights,
						use_max,
						black_init,
						
						exclusive_autosave,
						keep_Idata,
						real_log2,
						clone_models};
					
					if(master)
					{
						gd.addChoice("scanning pattern", pattern, scan_pattern);
						
						gd.addCheckboxGroup(4, 4, labels1, states1);
						gd.addNumericField("debug level", debug_level, 0);
						//gd.addMessage("0 .. none, 1 .. Pval.txt, 2 .. autosave, 3 .. Cachreport.txt, 4 .. job.bin");
						gd.addCheckbox("employ Smart_Atomizer", external_merit);
						gd.addCheckbox("lattice hopping", lattice_hopping);
					}
					gd.addCheckbox("update and display options.txt", display_options);
					gd.addChoice("Data input", titleArray, Datatitle);
					gd.addNumericField("Noisify (<=0 none)", noise_level,3);
					gd.addChoice("log Prob_map", ptableArray, PtableTitle);
					gd.addChoice("Beam profile", profileArray, beamProfileTitle);
					
					if(master)
					{
						gd.addChoice("Initial models", titleArrayPlus, Modeltitle);
						gd.addChoice("Action: ", Modii, Modus);
					}	
					else
					{	gd.addMessage("Action: " + Modus);}
					gd.addCheckbox("append another job (dataset)", another_job);
					if(!batch_run || (master && master_shows_dialog && (!suppress_further_dialog)) )
					{
						gd.showDialog();
						if (gd.wasCanceled())
						{	return;}
					}
					if(master)
					{
						timeoffset = System.currentTimeMillis();
						modelnum = (int) gd.getNextNumber();
						shadownum = (int) gd.getNextNumber();
						multiActivations= (int) gd.getNextNumber();
						multiActivations0 = multiActivations;
						
						modelactiveinit = (int) gd.getNextNumber();
						modelGraylevels = (int) gd.getNextNumber();
						
						if(modelGraylevels < 8)
						{
							IJ.error("ModelGraylevels is set to 8");
							modelGraylevels = 8;
						}

						modelBrightness = gd.getNextNumber();
						if(modelBrightness == 0.0)
						{
							IJ.log("brightness defaulted to -4");
							//for usual graphene the max carbon contrast is 4sigma
							modelBrightness = -4.0;
						}

						bondlength = (int) gd.getNextNumber();
						wobble = (int) gd.getNextNumber();
						initblocksize = (int) gd.getNextNumber();
						initblockdepth = (int) gd.getNextNumber();
						blocksize = initblocksize;
						if(blocksize == 0) {blocksize = 1;}
						if (blocksize < 0)
						{
							blocksize = -blocksize;
							flip_block = true;
						}
						new_blocksize = blocksize;
						skip_threads = (int) gd.getNextNumber();
						clone_weight = gd.getNextNumber();
					}
					if(( skip_threads > (max_threads-use_threads - 1) )) //can at most skip all but the last available workers
					{	skip_threads = 0;}
					if(master && (skip_threads > 0) )
					{
						if(!skipThreads(skip_threads))
						{	
							IJ.error("Error while skipping workers.");
							saveOptions(options, false);
							return;
						}
						skip_threads = 0;
					}
					initthreadnum = (int) gd.getNextNumber();
					threadnum = initthreadnum;
					if ( (threadnum < 1) || ( threadnum > (max_threads-use_threads) ) )
					{
						threadnum = max_threads - use_threads;
					}
					if(master)
					{
						thread_offset = (int) gd.getNextNumber();
						
						init_target_models = (int) gd.getNextNumber();
						target_models = init_target_models;
						init_target_wraps = (int) gd.getNextNumber();
						target_wraps = init_target_wraps;
					}
					init_use_frames = (int) gd.getNextNumber();
					use_frames = init_use_frames;
					if(master)
					{
						init_lattice_weight = gd.getNextNumber();
						initbenchmark = gd.getNextBoolean();
						if(batch_run)
						{	suppress_further_dialog = gd.getNextBoolean();}
					}
					shuffle_frames = gd.getNextBoolean();
					previous_seed = gd.getNextBoolean();
					if(!previous_seed)
					{
						seed = System.currentTimeMillis();
					}
					random = new Random(seed);
					if(master)
					{
						scan_pattern = gd.getNextChoice();
						alternating_scan_direction = gd.getNextBoolean();
						spiral = scan_pattern.equals("spiral");
						shuffle_update = scan_pattern.equals("noise");
						normalize_model_intensities = gd.getNextBoolean();
						blobs = gd.getNextBoolean();
						correlated_optimization = gd.getNextBoolean();
						quick_optimization  = gd.getNextBoolean();
						fixed_lattice_weight = gd.getNextBoolean();
						frozen_lattice = gd.getNextBoolean();
						use_weights = gd.getNextBoolean();
						next_use_weights = use_weights;
						fair_weights = gd.getNextBoolean();
						const_weights = gd.getNextBoolean();
						use_max = gd.getNextBoolean();
						next_use_max = use_max;
						black_init = gd.getNextBoolean();
						exclusive_autosave = gd.getNextBoolean();
						keep_Idata = gd.getNextBoolean();
						real_log2 = gd.getNextBoolean();
						clone_models = gd.getNextBoolean();
						debug_level = (int) gd.getNextNumber();
						external_merit = gd.getNextBoolean();
						lattice_hopping = gd.getNextBoolean();
					}
					Datatitle = gd.getNextChoice();
					noise_level = gd.getNextNumber();
					PtableTitle = gd.getNextChoice();
					beamProfileTitle = gd.getNextChoice();
					use_beamprofile = !beamProfileTitle.equals("<none>");
					
					if(use_beamprofile && (do_pix_EM || do_full_EM))
					{
						IJ.log("Error: beam proilfe " + beamProfileTitle + " cannot be used for any EM");
						return;
					}
					
					if( (do_pix_EM || do_full_EM) && ( (blocksize != 1) || (blockdepth != 1) || blobs) )
					{
						IJ.log("Error: EM can only be employed for 1x1x1 blocks withot blobs");
						return;
					}
					
					// the image we get the values from
					
					if(keep_Idata && duplicate_Idata)
					{
						//IJ.log("You should not keep AND duplicate Idata");
						duplicate_Idata = false;
					}
					
					display_options = gd.getNextBoolean();
					if(master)
					{
						if(external_merit)
						{
							sharedmerit = new SharedMerit();
							Smart_Atomizer.sharedmerit = this.sharedmerit;
							Smart_Atomizer.external_merit_source = getClass().getSimpleName();
							if(Smart_Atomizer.get_instance_count()==0)
							{
								IJ.doCommand("Smart_Atomizer");
							}
						}
						
						master_prefix = prefix;
						master_PtableTitle = PtableTitle;
						Modeltitle = gd.getNextChoice();
						//smooth_passes = (int) gd.getNextNumber();
						//Weighttitle = gd.getNextChoice();
						Modus = gd.getNextChoice();
						
						if(!options.isFile())
						{	
							saveOptions(options, false);
							IJ.log("created " + options.getPath() );
						}
						
						if( (!(PtableTitle.equals("<new>") || PtableTitle.equals("<none>") )) && ( Modeltitle.equals("avg. data") || Modeltitle.equals("none/noise") ) )
						{
							IJ.error("external Ptable does also require external models!");
							saveOptions(options, false);
							return;
						}
					}	
					
					if( Modus.equals(Modii[1]) || Modus.equals(Modii[2]))
					{	usenodes = new File( master_path + "usenodes.txt");}
					
					//only employ threads if we are not going to read them from usenodes.txt anyways
					if ( !(Modus.equals(Modii[2]) && usenodes.isFile() && usenodes.setReadable(true)) && !employThreads(threadnum) )
					{
						IJ.error("Error while employing workers.");
						saveOptions(options, false);
						return;
					}
					
					another_job = gd.getNextBoolean();
					if(display_options)
					{
						saveOptions(options, true);
						IJ.open(options.getPath());
						IJ.selectWindow(options.getName());
						saveOptions(options, false);
					}
					IJ.log("Action to perform: " + Modus);
					skip_stats = (!PtableTitle.equals("<new>") ) &&
								 (!( Modeltitle.equals("avg. data") || Modeltitle.equals("none/noise") ));
					/*
					if(!PtableTitle.equals("<new>") && (do_pix_EM || do_full_EM))
					{
						IJ.log("Error: external Ptable " + PtableTitle + " cannot be used for any EM");
						return;
					}
					*/
					if(master && Modus.equals(Modii[2]))
					{
						if(!loadJobrun(jobrun)) return;
					}
					
					if(another_job && (debug_level > 3) )
					{
						debug_level = 3;
						IJ.showMessage("Cannot write job.bin for parallell jobs\n" +
						"You could still create sepparate jobfiles and launch all at once");
					}
					
					minfornextm = 3*modelsize*modelsize / 4; //just in case we would use job files
					impWidth = modelsize; 
					impHeight = modelsize;
					if(master) //may be overriden by actual modelsize from Idata
					           //but needed if reading jobfiles
					{
						setHexcenter(modelsize/2, modelsize/2);
						minfornextm = 3*modelsize*modelsize / 4; //Number of real pixels
						hex_mask = new boolean[modelsize*modelsize];
						for(int ind = 0; ind < hex_mask.length ; ++ind)
						{	
							int q = ind%modelsize;
							int r = ind/modelsize;		
							hex_mask[ind] = iswithinHex(q,r,modelsize/2);		
						}
						
						
					}
					
					
					if(master && Modus.equals(Modii[2]))
					{
						job_stats_ok = readJobStats(jobStats);
						if(!job_stats_ok) 
						{
							IJ.log("failed to read " + jobStats.getPath() );
							return;
						}
					}
					
					if(!Modus.equals(Modii[2]))
					{		
						if(Datatitle.equals("<none>"))
						{
								IJ.log("Nothing to do without data");
								return;
						}
						if(noise_level > 0)
						{	//Job_Data will apply Poisson_Noise to Idata
							jobVec.get(job_id).noise_level = noise_level;
						}
						if(!Datatitle.equals("<generate>"))
						{	jobVec.get(job_id).Idata = WindowManager.getImage(Datatitle);}
						else
						{
							if(master)
							{
								datagen = new Data_Generator("Idata.tif", use_frames, bondlength);
								datagen.run( getClass().getSimpleName() );
								jobVec.get(job_id).Idata = WindowManager.getTempCurrentImage(); 
							}
							else
							{	jobVec.get(job_id).Idata = datagen.rerun("Idata", use_frames);}
							use_frames = jobVec.get(job_id).Idata.getStackSize();
							duplicate_Idata = false;
						}
						impWidth = jobVec.get(job_id).Idata.getWidth();   //get and check sizes
						impHeight = jobVec.get(job_id).Idata.getHeight();
						if(impHeight != impWidth)
						{
							IJ.error("Data Input must be square shape.");
							saveOptions(options, false);
							return;
						}
						modelsize = impWidth; //TODO only use modelsize hereafter
						modelactive = modelactiveinit;
						if(modelactive <= -60) //DUPLICATE CODE HERE
						{	modelactive = (-60*modelsize)/modelactive;}
						else if (modelactive < 0)
						{	modelactive = -modelactive * bondlength;}
						
						if( (modelactive > modelsize) || (modelactive < 2*bondlength) )
						{	modelactive = modelsize;}
						minfornextm = 3*modelsize*modelsize / 4; //Number of real pixels
						
					
						//normalize performance values (total is 1)
						IJ.log("Job" + job_id + " performance: " + normalize(jobVec.get(job_id).performances));
						jobVec.get(job_id).set_threadnum(threadnum); //creates new empty portions
					
					
						if(master)
						{
							setHexcenter(modelsize/2, modelsize/2);
							if(hex_mask == null || hex_mask.length!=modelsize*modelsize)
							{
								hex_mask = new boolean[modelsize*modelsize];
								for(int ind = 0; ind < hex_mask.length ; ++ind)
								{	
									int q = ind%modelsize;
									int r = ind/modelsize;		
									hex_mask[ind] = iswithinHex(q,r,modelsize/2);		
								}
							}
						}
						int actualstackSize = jobVec.get(job_id).Idata.getStackSize();
					
						stackSize = actualstackSize;
						if ( (use_frames > 0) && (use_frames < stackSize) )	
						{	stackSize = use_frames;}
						jobVec.get(job_id).set_use_frames(stackSize); //creates initial perm(utation)
						total_subframes += stackSize;
					
						//we use our own random generator for shuffling
						if (shuffle_frames)
						{
							shuffle(jobVec.get(job_id).perm);
							jobVec.get(job_id).cancel_portions(); //just to be safe
							//for(int i = 1; i < stackSize; ++i)
							//{	jobVec.get(job_id).swap_perm( i, random.nextInt(i) ); }
						}
						
						jobVec.get(job_id).init_stats(skip_stats); //need to know Idata and use_frames for that
														 //determines zipFrameLens
						datamin = jobVec.get(job_id).datamin;
						datamax = jobVec.get(job_id).datamax;
						datasum = jobVec.get(job_id).datasum;
						if(skip_stats)
						{
							datamean = 100;
							data_bg = 80;
							data_bg_std = 10;
							data_num_hexes = 10000;
							dataStdDev = 10;
						}
						else
						{
							datamean = jobVec.get(job_id).datamean;
							data_bg = jobVec.get(job_id).data_bg;
							data_bg_std = jobVec.get(job_id).data_bg_std;
							data_num_hexes = jobVec.get(job_id).data_num_hexes;
							dataStdDev = jobVec.get(job_id).dataStdDev;
							dataSumImg = jobVec.get(job_id).dataSumImg;
						}
						if(master)
						{	bigSumImg = new long[modelsize * modelsize];}
						if(!skip_stats)
						{
							for(int i = 0; i< bigSumImg.length; ++i)
							{	bigSumImg[i] += dataSumImg[i];}
						}
					
						if(!jobVec.get(job_id).fill_portions()) //assigns frames according to performances
						{	return;}							//calculates zipImgLen according to frame permutation
						
						IJ.log(Datatitle);
						
						{
							//Area per lattice hexagon 0.05238760873 (1.5 * 0.142 * 0.142 * Math.sqrt(3) )
							double actual_frame_area = 0.05238760873 * 0.25 * Math.pow( modelsize / bondlength, 2);
							IJ.log(stackSize + (shuffle_frames?" shuffled":"") +
								" frames, " + impWidth + "x" + impHeight +
								 " Pixels, " + actual_frame_area + " nm² each");
						}
						
						IJ.log("Data min, max, sum: " + datamin + ", " + datamax + ", " + datasum);
						
						
						if(skip_stats)
						{
							IJ.log("further stats are skipped!");
						}
						else
						{	
							IJ.log("Data heaxagons: " + data_num_hexes + "  background: " + data_bg + " stdDev: " + data_bg_std);
							IJ.log("Pixel mean +- stdDev: " + format.format(datamean) + " +- " + format.format(dataStdDev));
							{
								double num_pix = 0.75*Math.pow(modelsize,2);
								IJ.log("Frame mean +- stdDev:" + format.format(num_pix * datamean) + "+-" + format.format(dataStdDev * Math.pow(num_pix,0.5)) );
								//hexpixels per lattice hexagon         //Area per lattice hexagon (1.5 * 1.42 * 1.42 * Math.sqrt(3) )
								double dose = (3 * datamean * bondlength * bondlength) /  5.238760873 ;
								IJ.log("Dose Counts/A² : " + format.format(dose));
								double actual_area = 5.238760873 * stackSize * 0.25 * Math.pow( modelsize / bondlength, 2);
								IJ.log("Total area: " + actual_area + " A²   ( " + ( actual_area * Math.pow( 10, -8) ) + " ym² )" );
							}
						}
						modelmin = data_bg;
						modelmax = -1.0;
						
						if(modelBrightness < 0)
						{	
							double brightness_unit = (datamean-modelmin);
							modelmin -= brightness_offset*brightness_unit; //barely sufficent for 8 ring
							brightness_unit = (datamean-modelmin);
							modelmax = modelmin - (modelBrightness * brightness_unit);
							//modelmax = modelmin - (modelBrightness * (1.0-brightness_offset)) * (datamean-modelmin);
							//modelmin += (modelBrightness * brightness_offset) * (datamean - modelmin );
							if(modelmin < 0)
							{
								modelmax -= modelmin;
								modelmin = 0;
							}	
						}
						else
						{	modelmax = modelmin + modelBrightness;}
						//modelmin = (double)datamin; //now we have background for that
						///DEBUG TESTING EM
						///modelmin = 635.42;
						///modelmax = 2753.85;
						
						IJ.log("Model min, max, levels : " + format.format(modelmin) + ", " + format.format(modelmax) + ", " + modelGraylevels);
						jobVec.get(job_id).modelmin = modelmin;
						jobVec.get(job_id).modelmax = modelmax;
						
						if(!(PtableTitle.equals("<new>") || PtableTitle.equals("<none>")))
						{
							int new_width = (1+datamax - datamin);
							pvaloffset = 0.0;
							Ptable = WindowManager.getImage(PtableTitle);
							if( ( Ptable.getWidth() == new_width ) &&
							    (Ptable.getHeight() == (modelGraylevels) ) )
							{	
								jobVec.get(job_id).Ptable = Ptable;
								IJ.log("Ptable: " + PtableTitle);	
							}
							else
							{
								
								if( ( Ptable.getWidth() >= new_width ) &&
									(Ptable.getHeight() == modelGraylevels ) )
								{
									IJ.log("cropping " + PtableTitle);
									ImageStack bigPtableSt = Ptable.getStack();
									ImageStack newPtableSt = bigPtableSt.crop(datamin,0,0,new_width,modelGraylevels,1);
									Ptable.setStack(newPtableSt);
									Thread.sleep(10);
									Ptable.updateAndRepaintWindow();
									Thread.sleep(10);
									Ptable.changes = false; //suppress save changes Dialog
									jobVec.get(job_id).Ptable = Ptable;
								}
								else
								{
									
									if(skip_stats)
									{
										saveOptions(options,true);
										IJ.error("Bad Luck: stats were skipped and Ptable could not be cropped, either!");
										return;
									}
									IJ.error("Dimensions of Ptable dont match!, defaulting to Poisson PDF");
									PtableTitle = "<new>";
								}
							}
						}
						
						if(PtableTitle.equals("<new>"))
						{
							jobVec.get(job_id).autopvaloffset(pvaloffset, autoptableoffset);
							//create lookup table for probabilites
							jobVec.get(job_id).createPtable(modelGraylevels);
							PtableTitle = jobVec.get(job_id).Ptable.getTitle();
						}
						
						if( use_beamprofile )
						{
							beamProfile = WindowManager.getImage(beamProfileTitle);
							beamsize = beamProfile.getWidth();
							beamarea = beamsize * beamsize;
							if(beamsize != beamProfile.getHeight())
							{
								IJ.log("WARNING: cannot use non squared " + beamProfileTitle);
								use_beamprofile = false;	
							}
							if(beamProfile.getType() != ImagePlus.GRAY16)
							{
								IJ.log("WARNING: " + beamProfileTitle + " is not GRAY16 formated");
								use_beamprofile = false;
							}
						}
						
						
					}
					
					
					//////////////
					// model greyval (int) = modelGraylevels * ((actual value) - modelmin)/(modelmax-modelmin)
					//
					// actual value = modelmin + (modelmax-modelmin)*( (float)(model grey val) / (float)modelGraylevels )

					if(master)
					{
						m_Datatitle = Datatitle;
						m_initthreadnum = initthreadnum;
						m_init_use_frames = init_use_frames;
						m_previous_seed = previous_seed;
						m_sourcepath = sourcepath;
						m_seed = seed;
						///TODO is that check really needed here!?
						if( (modelactive > modelsize) || (modelactive < 2 * bondlength) )
						{	modelactive = modelsize;}
						IJ.log("blobs (very expensive!): " + blobs);
						IJ.log("Scanning pattern: " + scan_pattern);
						IJ.log("Models: " + Modeltitle);
						if((modelnum < 1) && ( Modeltitle.equals("none/noise") || Modeltitle.equals("avg. data")) )
						{
							IJ.error("Cannot create " + modelnum + " defaulting to " + defaultmodelnum );
							modelnum = defaultmodelnum;
						}
						if ( Modeltitle.equals("none/noise") )
						{
							Modeltitle = "NewModels.tif";
							new_models = true;
							Models = NewImage.createShortImage( Modeltitle , modelsize , modelsize , modelnum, NewImage.FILL_BLACK);
						}
						else if ( Modeltitle.equals("avg. data") ) 
						{
							if(Modus.equals(Modii[2]) && !job_stats_ok) 
							{	
								IJ.log("cannot re-create avg data, failed to open JobStats.txt");
								return;
							}
							Modeltitle = "AvgData.tif";
							Models = NewImage.createShortImage( Modeltitle , modelsize , modelsize , modelnum, NewImage.FILL_BLACK);
							ImageStack localModelSt = Models.getStack();
							double inv_Stsize = Modus.equals(Modii[2]) ? (double)1.0/dataStSize0 : (double)1.0/(total_subframes);
							
							//initialize first model with scaled average
							if(!Modus.equals(Modii[2]))
							{
								data_bg0  = jobVec.get(default_job).data_bg;
								data_bg_std0  = jobVec.get(default_job).data_bg_std;
								modelmin0 = jobVec.get(default_job).modelmin;
								modelmax0 = jobVec.get(default_job).modelmax;
							}
							
							{
								
								//float[] probs = null;
								short[] pixels = (short[]) localModelSt.getPixels(1);
								//if( jobVec.get(default_job).Ptable != null)
								//{	probs  = (float[]) jobVec.get(default_job).Ptable.getStack().getPixels(1);}
								int mean = 0;
								for(int i = 0; i < pixels.length; ++i)
								{	
									short val = (short) Math.round( modelGraylevels *
									( (double)bigSumImg[i] * inv_Stsize - data_bg0 - modelmin0 ) / (modelmax0 - modelmin0) );
									if( val >= modelGraylevels) {val = (short) (modelGraylevels - 1);}
									if( val < 0) {	val = 0;}
									pixels[i] = val;
									mean += val;
								}
								mean = (mean + pixels.length/2) / pixels.length;
								
								/*for(int j = 0; j < smooth_passes; ++j)
								{
									for(int i = 0; i < pixels.length; ++i)
									{	
										int q = i % modelsize;
										int r = i / modelsize;
										if( iswithinHex(q,r,modelsize/2))
										{	pixels[i] = (short)( (2 + 4*pixels[i] + (short)mean)/5 );}
									}
								}*/
								//in most cases only the first one is needed
								//but this one also works if Weights.txt already exists 
								for(int m = 2; m <= modelnum; ++m)
								{	
									short[] pix = (short[]) localModelSt.getPixels(m);
									for(int i = 0; i < pixels.length; ++i)
									{	pix[i] = pixels[i]; }
								}
								Models.setStack(localModelSt);
							}
						}
						else
						{
							Models = WindowManager.getImage(Modeltitle);
							if(Models == null)
							{
								IJ.log("could not find " + Modeltitle);
								saveOptions(options, false);
								return;
							}
							if(Models.getType() != ImagePlus.GRAY16 || modelsize != Models.getWidth() || modelsize != Models.getHeight()  )
							{
								IJ.error( Modeltitle + " must be a stack of GRAY16 " + modelsize +"x"+modelsize+" images ");
								saveOptions(options, false);
								return;
							}
						}
						if(modelnum < 1) { modelnum = Models.getStackSize();}
						
						nextModel = NewImage.createShortImage( "nextModel" , modelsize , modelsize , 1, NewImage.FILL_BLACK);
						{
							short[] nextpixels = (short[])nextModel.getStack().getPixels(1);
							short[] pixels = (short[])Models.getStack().getPixels(1);
							for(int i = 0; i < pixels.length; ++i)
							{
								nextpixels[i] = pixels[i];
							}
						}
						if(sym_scaling)
						{
							nextMatcher = new Diff_Matcher(nextModel, 1, 1, bondlength, (median_filtering?modelGraylevels:0) );
							//nextMatcher.exponent = sym_fac_exponent;
						}
						//initialize models only AFTER weigts have been initialized
						if( (target_models < 0) || (target_models > modelnum) ) {target_models = modelnum;}

						if(multiActivations > modelnum)
						{	multiActivations = modelnum;}
						if(multiActivations < 0)
						{
							multiActivations = modelnum/(-multiActivations);
						}
						if(multiActivations == 0)
						{	multiActivations = modelnum;}


						if(init_lattice_weight < 0.0)
						{	init_lattice_weight = -init_lattice_weight/modelnum;}
						else if(init_lattice_weight > 1.0)
						{	init_lattice_weight -= Math.floor(init_lattice_weight);}
						
						
						if(blocksize > modelsize/2)
						{
							blocksize = modelsize/2;
						}
						// <= 0 seems mor practical to preseve model0 as empty graphene
						reverse_block_depth = (initblockdepth <= 0);
						
						if(initblockdepth < 0)
						{
							blockdepth = -initblockdepth;
						}
						else if( (initblockdepth == 0) || (initblockdepth > modelnum) )
						{	
							blockdepth = modelnum;
						}
						else
						{	blockdepth = initblockdepth;}

						if( (blockdepth > modelnum) || (blockdepth < 1))
						{
							blockdepth = modelnum;
						}
						new_blockdepth = blockdepth;
						blockarea = blocksize * blocksize;
						blockvolume = blockarea * blockdepth;
						IJ.log("Blocksize: " + blocksize + " Blockdepth: " + blockdepth);

						{
							int num_blocks = (modelsize * modelsize) / (blocksize * blocksize) ;
							//int seed = random.nextInt(num_blocks);
							//There may be a bug from https://gcc.gnu.org/bugzilla/show_bug.cgi?id=8481
							//IT SHOULD HAVE BEEN FIXED IN 2002 !!!
							IJ.log("Number of Blocks: " + num_blocks);
							/*
							while( (seed < 0) && (seed >= num_blocks) && (num_blocks > 0))
							{
								seed = random.nextInt(num_blocks);
								IJ.log("invalid seed, new seed: " + seed);
							}
							*/
							permutation = new int[modelnum][num_blocks];
							init_perm(permutation);
							//permutator = new Permutator(num_blocks,random.nextInt(num_blocks));
						}
						//array to count unsuccessfull attempts at pixel modification in each model image
						unsc = new int[modelnum];
						dead = new ModelState[modelnum];

						dead[0] = ModelState.LATTICE;	
						for (int i = 1; i < modelnum; ++i)
						{	dead[i] = ModelState.RESET; }//Will be initialized by default
						
						
						IJ.log("Weights: " + Weighttitle);
						IJ.log("use weights: " + use_weights + ( (use_weights && const_weights)?" (const)":"") + ( (use_weights && aggressive_weights)?" (aggressive)":"") );
						IJ.log("enable pixel EM: " + do_pix_EM);
						IJ.log("enforce pixel EM: " + force_pix_EM);
						IJ.log("full Model EM: " + do_full_EM);
						IJ.log("collect Histogramms: " + collect_Histos);
						
						modelweight = new double[modelnum];
						modelweightupd = new double[modelnum];
						sym_factors = new double[modelnum];
						cases = new int[modelnum];
						if ( !initialize_weights(master_path, Weighttitle) )
						{	return;}
						//dont save them before models have been initalized
						{
							File peek = new File(master_path + "Diff.tif");
							if(peek.isFile())
							{
								Diffimg = new ImagePlus(peek.getPath());
								IJ.log ("opened Diff.tif");
							}
							if(Diffimg == null || Diffimg.getStackSize()!=modelnum || Diffimg.getWidth() != modelsize || Diffimg.getHeight() != modelsize)
							{
								Diffimg = NewImage.createShortImage("Diff.tif", modelsize , modelsize , modelnum, NewImage.FILL_BLACK);
							}
							DiffSt = Diffimg.getStack();

							if(! peek.isFile())
							{
								for (int v=1; v<= DiffSt.getSize(); ++v)
								{
									short[] diffpx = (short[])DiffSt.getPixels(v);
									for (int i=0; i<diffpx.length; ++i)
									{
										diffpx[i]=(short)32768;
									}
								}
							}
							Diffimg.show();
						}
						IJ.log(Modeltitle + " (WxHxD) " + modelsize + "x" + modelsize + "x" + modelnum);
						IJ.log("Model active region:" + modelactive);
						//initialize Models AFTER weights and Diff exist
						if(!Modus.equals(Modii[2])) //otherwise already read from jobStats.txt
						{
							data_bg0  = jobVec.get(default_job).data_bg;
							datamean0 = jobVec.get(default_job).datamean;
							modelmin0 = jobVec.get(default_job).modelmin;
							modelmax0 = jobVec.get(default_job).modelmax;
							
							//Funnily enough we have to pick one job, but it doesnot matter which one.
						}	
						
						
						modelmean0 = ( ( (double)modelGraylevels) * ( datamean0 - data_bg0 - modelmin0 ) / ( modelmax0 - modelmin0 ) );
						initialize_models( Models, Diffimg, modelGraylevels, modelmean0);
						ModelsSt = Models.getStack();
						Models.show();
						//Models.updateAndRepaintWindow();
						Diffimg.show();
						//Diffimg.updateAndRepaintWindow();
						//save and normalize weights AFTER models have been initialized
						normalizeWeights(modelweight); //will crash if fixed lattice weight != 1
						{
							
							File peek = new File(master_path + "Weights.txt");
							if(!peek.isFile())
							{
								saveWeights(master_path);
							}
						}
						
						wrapnum = 0;
						histwidth = modelnum * ( 1 + modelsize);
						histheight = 3 * modelsize + 2;
						
						{
							boolean loading_failed = true;
							File peek = new File(master_path + "Model_history.tif");
							if(!peek.isFile())
							{	peek =  new File(master_path + "Model history.tif");} //old_version
							
							if ( peek.isFile() )
							{
								ImagePlus loadedHistimg = new ImagePlus(peek.getPath());
								wrapnum = loadedHistimg.getStackSize() - 1;
								if(loadedHistimg.getWidth() == histwidth )
								{
									if (loadedHistimg.getHeight() != histheight)
									{
										IJ.log("upgrading Model_history.tif from 2 rows to 3 rows");
										Histimg = NewImage.createFloatImage("Model_history.tif",histwidth, histheight,wrapnum+1,NewImage.FILL_BLACK);
										int oldHeight = loadedHistimg.getHeight();
										for(int w = 1; w <= wrapnum + 1; ++w)
										{
											float[] oldpix = (float[])loadedHistimg.getStack().getPixels(w);
											float[] newpix = (float[])Histimg.getStack().getPixels(w);
											for(int oldi=0; oldi < oldpix.length; ++oldi)
											{
												int oldy = oldi / histwidth;
												int newy = oldy;
												if((oldy > modelsize) && (oldy < 2 * modelsize + 1) ) //we are in the second row
												{	newy = oldy + modelsize + 1;} //move down one row
												int newi = newy * histwidth + (oldi % histwidth);
												newpix[newi] = oldpix[oldi];
											}	
										}
										float[] histpix = (float[])Histimg.getStack().getPixels(wrapnum+1);
										for(int m = 0; m < modelnum; ++m)
										{
											short[] diffpix = (short[])DiffSt.getPixels(m + 1);
											for(int y = 0; y < modelsize; ++y)
											{
												for(int x = 0; x < modelsize; ++x)
												{
													int intdiff = ( diffpix[x+modelsize*y] & 0xffff );
													double diffval = ( (double)(intdiff - 32768) ) * 
																	 (modelmax0-modelmin0) / 
																	 ( (double)modelGraylevels ); //actual intensity
													histpix[x+m*(modelsize+1) + histwidth*(y+modelsize+1)] = (float)diffval;
												}
											}
										}
									}
									else
									{
										IJ.log ("opened Model_history.tif wrapnum: " + wrapnum);
										Histimg = loadedHistimg;
									}
									loading_failed = false;
									append_Pval = true; //reset Pval unless ModelHistory is present
									if(target_wraps > 0)
									{	target_wraps += wrapnum;}
								}
								else
								{
									IJ.log("Dimensions dont match, reinitializing Model_history.tif");
									if(target_wraps > 0)
									{	target_wraps += wrapnum;}
								}
							}

							if(loading_failed)
							{
								Histimg = NewImage.createFloatImage("Model_history.tif",histwidth, histheight,1,NewImage.FILL_BLACK);
								ImageStack HistimgSt = Histimg.getStack();
								float[] histpix = (float[]) HistimgSt.getPixels(1);
								int modelarea = modelsize * modelsize;
								for (int m = 0; m < modelnum; ++m)
								{
									short[] pixels = (short[])ModelsSt.getPixels(m+1);
									short[] diffpix = (short[])DiffSt.getPixels(m+1);
									for(int ind = 0; ind < modelarea; ++ind)
									{
										int x = ind % modelsize;
										int y = ind / modelsize;
										float val = /*(float)modelmin0+*/((float)pixels[x+modelsize*y])*
											((float)(modelmax0-modelmin0))/((float)modelGraylevels);
											
										int intdiff = ( diffpix[x+modelsize*y] & 0xffff );
										double diffval = ( (double)(intdiff - 32768) ) * 
														(modelmax0-modelmin0) / 
														( (double)modelGraylevels ); //actual intensity
										
										histpix[x+m*(modelsize+1)+histwidth*y] = val;
										histpix[x+m*(modelsize+1)+histwidth*(y+modelsize+1)] = (float)diffval;
										histpix[x+m*(modelsize+1)+histwidth*(y+2*modelsize+2)] = (float)diffval;										
									}
								}
							
							}
						}
						//ImageProcessor Histimgp = Histimg.getProcessor();
						Histimg.show();
						//Histimg.updateAndRepaintWindow();	
						File peek = new File(master_path + "Updates.tif");
						boolean updates_loaded = false;
						if (peek.isFile())
						{
							ImagePlus loadedUpdates = new ImagePlus(peek.getPath());
							if(	loadedUpdates.getWidth() == modelsize &&
								loadedUpdates.getHeight() == modelsize &&
								loadedUpdates.getStackSize() == modelnum)
							{
								Updates = loadedUpdates;
								updates_loaded = true;
								IJ.log("loaded " + loadedUpdates.getTitle() );
							}
							else
							{
								IJ.log("WARNING wrong dimensions in " + loadedUpdates.getTitle() );
							}	
						}
						
						if(!updates_loaded)
						{	
							Updates = NewImage.createShortImage( "Updates.tif" , modelsize , modelsize , modelnum, NewImage.FILL_BLACK);
							IJ.log("created new " + Updates.getTitle());	
						}
						Updates.show();
						UpdatesSt = Updates.getStack();
						update_pixels = (short[])UpdatesSt.getPixels(1);
					
						peek = new File(master_path + "Activity.tif");
						boolean activity_loaded = false;
						if (peek.isFile())
						{
							ImagePlus loadedActivity = new ImagePlus(peek.getPath());
							if(	loadedActivity.getWidth() == modelsize &&
								loadedActivity.getHeight() == modelsize &&
								loadedActivity.getStackSize() == modelnum)
							{
								Activity = loadedActivity;
								activity_loaded = true;
								IJ.log("loaded " + loadedActivity.getTitle() );
							}
							else
							{
								IJ.log("WARNING wrong dimensions in " + loadedActivity.getTitle() );
							}	
						}
						if(!activity_loaded)
						{	
							Activity = NewImage.createShortImage( "Activity.tif" , modelsize , modelsize , modelnum, NewImage.FILL_BLACK);
							IJ.log("created new " + Activity.getTitle() );
						}
						
						
						ActivitySt = Activity.getStack();
						if(!activity_loaded || reset_activity)
						{
							int active_radius = (wrapnum==0) ? modelactive/2 : modelsize/2;
							for(int m = 0; m < modelnum ; ++m)
							{
								int mm = (modelnum - m - 1)/multiActivations;
								{ mm = (mm % 14)+1 ;} //min 1 max 14 
								
								//System.out.println("model: " + m + " launch_delay: " + mm);
								
								activity_pixels = (short[])ActivitySt.getPixels(m+1);
								for(int i = 0; i < activity_pixels.length ; ++i)
								{	
									activity_pixels[i] = (short)32768; 
									int q = i % modelsize;
									int r = i / modelsize;
									int mmm = mm;
									if(iswithinHex(q,r,modelactive/2))
									{	
										--mmm;
										if(mmm==0)
										{	activity_pixels[i] += 66;}	
									} 
									if(hex_mask[i])//(iswithinHex(q,r,modelsize/2))
									{	activity_pixels[i] -= 2*Math.pow(2,mmm);}
								}
							}
						}
						Activity.show();
						//Activity.updateAndDraw();
						if( Modus.equals(Modii[1]) )
						{
							saveAllImgs( master_path);
							//we dont have yet any Pvalues that would need a backup
						}
						if( Modus.equals(Modii[3]) )
						{
							modelManipulation(true);
							saveAllImgs( master_path);
							saveWeights(master_path);
							if(IJ.showMessageWithCancel(getClass().getSimpleName()+ " asks",
                                            "proceed with reconstruction?") )
                            {                
								Modus=Modii[0];
								IJ.log("running EXPERIMENTAL FEATURE");
								IJ.log("next Action: " + Modus);
								histwidth = Histimg.getWidth();
							}
							else
							{
								return;
							}
						}
					
						
					}
					
					date.setTime(System.currentTimeMillis());
					if(master)
					{
						if( Modus.equals(Modii[1]) )
						{
							usenodes_stream = new java.io.FileOutputStream(usenodes.getPath(), false );
							usenodes_writer = new java.io.PrintWriter(usenodes_stream);
							usenodes_writer.println("nodes required for these job files on " + dateFormat.format(date));
							usenodes_writer.flush();
							sendjobstream = new java.io.FileOutputStream(sendjob.getPath(), false );
							sendjobwriter = new java.io.PrintWriter(sendjobstream);
							sendjobwriter.println("#!/bin/bash");
							sendjobwriter.println("#send jobfiles to all workers");
							sendjobwriter.println("#" + dateFormat.format(date) );
							sendjobwriter.flush();
						}
						else
						{
							nukestream = new java.io.FileOutputStream(nuke.getPath(), false );
							nukewriter = new java.io.PrintWriter(nukestream);
							nukewriter.println("#!/bin/bash");
							nukewriter.println("#nuke all workers of " + user);
							nukewriter.println("#VERIFY that these are YOUR zombie WORKERS");
							nukewriter.println("#" + dateFormat.format(date) );
							nukewriter.flush();
							
						}
					}
					//launch worker threads
					if(Modus.equals(Modii[2]))
					{
						if(usenodes.isFile() && usenodes.setReadable(true) )
						{	
							usenodes_sc = new Scanner(new FileReader(usenodes));
							usenodes_sc.nextLine();
							if(check_free_nodes)
							{	IJ.log("usenodes.txt will override worker availability");}
							threadVec.clear();
							threadadd = 0;
							while(usenodes_sc.hasNextLine())
							{	
								String line = usenodes_sc.nextLine();
								threadVec.add( new ThreadData());
								threadVec.get(threadadd).command = line;
								Scanner lsc = new Scanner(line);
								int start = line.lastIndexOf("job") + 3;
								int end = line.lastIndexOf(".bin");
								int worker = Integer.parseInt( line.substring(start,end) );
								threadVec.get(threadadd).worker = worker;
								lsc.close();
								++threadadd;
							}
							firstthread = 0;
							use_threads = threadadd;
							max_threads = threadadd;
							if( max_threads > p.length ) // create new arrays
							{	//we need larger arrays for more workers
								
								IJ.log("ERROR: you are exceeding the limit of " + MAX_THREADS + " workers");
								return;
								//THIS does not work. Probably streams or something cannot be copied use hard limits instead
								/*
								java.lang.Process p2[] = new java.lang.Process[max_threads];
								java.io.InputStream is2[] = new java.io.InputStream[max_threads];
								java.io.BufferedReader reader2[] = new java.io.BufferedReader[max_threads];
								java.io.OutputStream os2[] = new java.io.OutputStream[max_threads + 1];
								java.io.PrintWriter pwos2[] = new java.io.PrintWriter[max_threads + 1];
								boolean can_write2[] = new boolean[max_threads + 1];
								// copy old values to beginning of new array
								for(int i = 0; i < p.length; ++i ) 
								{
									p2[i] = p[i];
									is2[i] = is[i];
									reader2[i] = reader[i];
									os2[i] = os[i];
									pwos2[i] = pwos[i];
									can_write2[i] = can_write[i];
								}
								//set pointers to new arrays, JAVA VM will garbage collect the old arrays
								p = p2;
								is = is2;
								reader = reader2;
								os = os2;
								pwos = pwos2;
								can_write = can_write2;
								*/ 
							}
						}
						IJ.log("broadcasting files for job" + job_id);
						Process ps2 = rt.exec( sendjob.getPath() );
						ps2.waitFor();
						ps2.destroy(); 
					}
					
					//launch or reserve workers
					for (int i = firstthread; i < use_threads; ++i)
					{
						int worker = threadVec.get(i).worker + thread_offset;
						//JAJA bad Performance TODO that test outside the for loop
						if(Modus.equals(Modii[0]) || (Modus.equals(Modii[2]) && usenodes.isFile() && usenodes.setReadable(true)) ) 
						{
							p[i] = rt.exec(threadVec.get(i).command);  //launch thread i
							os[i] = p[i].getOutputStream(); //get process i's input stream to os[i]
							is[i] = p[i].getInputStream();  //get process i's output to is[i]
							reader[i] = new java.io.BufferedReader(new InputStreamReader(is[i]));
						}
						else if(Modus.equals(Modii[2]))
						{
							if(usenodes_sc == null)
							{
								if(threadVec.get(i).command.startsWith("/home"))
								{
									p[i] = rt.exec( threadVec.get(i).command + " /home/" + user + "/job" + worker + ".bin" );
								}
								else
								{
									p[i] = rt.exec(threadVec.get(i).command + " job" + worker + ".bin");
								}
							}
							else
							{
								System.out.println("Error: this code should be dead");
								p[i] = rt.exec(threadVec.get(i).command);
							}
							
							os[i] = p[i].getOutputStream(); //get process i's input stream to os[i]
							is[i] = p[i].getInputStream();  //get process i's output to is[i]
							reader[i] = new java.io.BufferedReader(new InputStreamReader(is[i]));
						}
						else if(Modus.equals(Modii[3]))
						{
							//dont launch any workes
							int load = (int) threadVec.get(i).load;
							int perf = (int) (threadVec.get(i).performance/threadVec.get(i).load);
							if(threadVec.get(i).command.startsWith("/home"))
							{
								usenodes_writer.println(threadVec.get(i).command + " /home/" + user + "/job" + worker + ".bin");
							}
							else
							{
								usenodes_writer.println(threadVec.get(i).command + " job" + worker + ".bin");
							}
							usenodes_writer.flush();
							
							IJ.log("creating " + parentpath + "job" + worker + ".bin");
							os[i] = new java.io.DataOutputStream(new java.io.FileOutputStream(parentpath + "job" + worker + ".bin", false));
							if( threadVec.get(i).command.startsWith("/home") )
							{
								//System.out.println("LOCAL:" + threadVec.get(i).command);
								sendjobwriter.println("cp " + parentpath + "job" + worker + ".bin /home/" + user + "/job" + worker + ".bin" );
								sendjobwriter.flush();	
							}
							else
							{
								Scanner sc = new Scanner(threadVec.get(i).command);
								//System.out.println("REMOTE:" + threadVec.get(i).command);
								sc.next();
								sc.next();
								int port = sc.nextInt();
								String host = sc.next();
								sendjobwriter.println("scp -P " + port + " " + parentpath + "job" + worker + ".bin " + host + ":");
								sendjobwriter.flush();
							}
						}
						else
						{
							IJ.error("Undefined Action: " + Modus);
							return;
						}
						pwos[i] = new java.io.PrintWriter(os[i]);
						can_write[i] = threadVec.get(i).can_write;
					}
					streamnum = use_threads;
					if(debug_level > 3)
					{
						IJ.log("writing " + master_path + "job.bin");
						java.io.FileOutputStream jobstream = new java.io.FileOutputStream(master_path + "job.bin", false);
						os[use_threads] =  new java.io.DataOutputStream(jobstream);
						pwos[use_threads] = new java.io.PrintWriter(jobstream);
						can_write[use_threads] = false; //dont requ
						++streamnum;
					}
					
					if(master && (debug_level > 0) )
					{
						File peek = new File( master_path + "Pval.txt");
						if(peek.isFile() && append_Pval)
						{
							IJ.log("appending " + master_path + "Pval.txt");
						}
						else
						{
							IJ.log("creating/overwriting " + master_path + "Pval.txt");
							looptotal = 0;
						}	
						Pvalstream = new java.io.FileOutputStream(master_path + "Pval.txt", append_Pval); //append new values
						Pvalhistory = new java.io.PrintWriter(Pvalstream);
					}
					
					
					date.setTime(System.currentTimeMillis());
					if(!Modus.equals(Modii[2]))
					{
						//begin with stuff that is sent to all worker threads
						for (int i=firstthread; i < streamnum; ++i)
						{
							//send all relevant variables
							
							int job = 0;
							int rank = 0;
							int worker = 0;
							if(i < use_threads)
							{
								job = threadVec.get(i).job;
								rank = threadVec.get(i).rank;
								worker = threadVec.get(i).worker + thread_offset;
							}
							
							pwos[i].println("ThreadNum(" + worker + ")");
							pwos[i].println("#" + jobVec.get(job).Idata.getTitle());
							pwos[i].println("#" + dateFormat.format(date));
							pwos[i].println("FrameCount(" + jobVec.get(job).portions[rank] + ")");
							pwos[i].println("ZipImgLen(" + jobVec.get(job).zipImgLen[rank]+ ")");
							pwos[i].println("PvalOffSet(" + jobVec.get(job).pvaloffset + ")");
							pwos[i].println("ModelSize(" + modelsize + ")");
							pwos[i].println("ModelNum(" + modelnum + ")");
							pwos[i].println("ShadowNum(" + shadownum + ")");
							pwos[i].println("ModelVal(" + modelGraylevels + ")");
							pwos[i].println("DataMax(" + jobVec.get(job).datamax + ")");
							pwos[i].println("DataMin(" + jobVec.get(job).datamin + ")");
							pwos[i].println("BondLength(" + bondlength + ")");
							pwos[i].println("Translations(" + (lattice_hopping?1:0) +")");
							pwos[i].println("Wobble(" + wobble + ")");
							pwos[i].println("debug_level(0)");
							pwos[i].println("pvalscaling(" + pvalscaling + ")");
							pwos[i].println("use_weights(" + (use_weights?1:0)  + ")");
							pwos[i].println("use_max(" + (use_max?1:0) +")");
							pwos[i].println("black_init(" + (black_init?1:0)  + ")");
							pwos[i].println("real_log2(" + (real_log2?1:0)  + ")");
							pwos[i].println("use_profile(" + (use_beamprofile?1:0) + ")");
							pwos[i].println("BeamSize(" + beamsize + ")");
							
							if(i==use_threads)
							{
								pwos[i].println("Reporting(1)");
							}
							else
							{
								pwos[i].println("Reporting(0)");
							}
							pwos[i].println("Run(1)");
							pwos[i].flush();
						}
						//IJ.log("sent out all initalizers: " + use_threads);
						//forward status/error reports of workers to log
					}
				
					s = "";	
					if( !Modus.equals(Modii[2]) )
					{
					    //messages like
					    //Worker0: FLP SUM(768) with cc_table (black init)  ( < 1750 MB) with pid: 3244 on cpu: 0 at Linux, weights: 1
						for (int i=firstthread; i<use_threads; ++i)
						{
							//workers report either from reading stdin or the jobfile
							boolean retry = false;
							do 
							{
								s = reader[i].readLine();
								if(s == null)
								{	
									System.out.println("no reply from worker" + i + " aka " + threadVec.get(i).command);
									IJ.error("no reply from worker" + i + " aka " + threadVec.get(i).command +
									"\n Click OK to try a second time");
									if(!retry)
									{
										retry = true;
										continue;
									}
									else
									{	return;}
								}
								
								if(!s.startsWith("Worker")) 
								{
									System.out.println(threadVec.get(i).command + " " +
									"job: " + threadVec.get(i).job +
									" worker: " + i + " : " + s);
								}
							} 
							while(!s.startsWith("Worker") && reader[i].ready() );
							IJ.log(threadVec.get(i).command);
							IJ.log(s);
							
							Scanner sc = new Scanner(s);
							while ( sc.hasNext() && !sc.next().equals("pid:") ) {};
							int pid = sc.hasNext()?sc.nextInt():-42;
							String s2 = threadVec.get(i).command;
							
							if(s2.startsWith("/home"))
							{
								nukewriter.println("kill -9 -pid " + pid);
							}
							else
							{
								//first four words "ssh -p 22xx cluster "
								Scanner lsc = new Scanner(s2);
								for(int j = 0; j < 4; ++j)
								{	nukewriter.print(lsc.next() + " ");}
								lsc.close();
								nukewriter.println(" kill -9 -pid " + pid);
							}
						}
					}
					
					if(!Modus.equals(Modii[2]))
					{
						if(use_beamprofile)
						{
							for (int i=firstthread; i<streamnum; ++i)
							{
								int job = 0;
								if( i < use_threads)
								{
									job = threadVec.get(i).job;
								}
								
								
								///THIS DOES NOT PREVENT CLUSTERHEAD FROM CRASHING DURING SUBMISSIONS
								/*if( (i != firstthread) && (i != 0) && !Modus.equals(Modii[1]) )
								{	
									System.out.println("sleep 100ms before sending data to " + threadVec.get(i).command);
									Thread.sleep(100);
								}*/
								pwos[i].println("updateBeamProfile()");
								pwos[i].flush();
								pwos[i].println("BeginBinary(" + (2 * beamarea) + ")");
								pwos[i].flush();
							
								int checksum = 0;
								ImageStack beamProfileSt = beamProfile.getStack();
								short[] beampix = (short[])beamProfileSt.getPixels(1);
								for(int ind = 0; ind < beamarea; ++ind)
								{
									short val = (short) (beampix[ind]&0xFFFF-32768);
									checksum += val;
									/*
									System.out.print("" + val);
									if((ind+1)%beamsize == 0)
									{
										System.out.print("\n");
									}
									else
									{
										System.out.print("\t");
									}
									*/ 
									send16bitint(os[i],val);
								}
								os[i].flush();
								pwos[i].println("EndBinary(" + checksum + ")");
								pwos[i].flush();
							}
						}
						
						
						
						
						
						for (int i=firstthread; i<streamnum; ++i)
						{
							int job = 0;
							if( i < use_threads)
							{
								job = threadVec.get(i).job;
							}
							datamax = jobVec.get(job).datamax;
							datamin = jobVec.get(job).datamin;
							
							///THIS DOES NOT PREVENT CLUSTERHEAD FROM CRASHING DURING SUBMISSIONS
							/*if( (i != firstthread) && (i != 0) && !Modus.equals(Modii[1]) )
							{	
								System.out.println("sleep 100ms before sending data to " + threadVec.get(i).command);
								Thread.sleep(100);
							}*/
							pwos[i].println("updatePtable()");
							pwos[i].flush();
							pwos[i].println("BeginBinary(" + (8 * (1+datamax-datamin) * modelGraylevels) + ")");
							pwos[i].flush();
						
							double fchecksum = 0.0;
							int threadsum = 0;
							Ptablep = jobVec.get(job).Ptable.getProcessor();
							for (int datav=0; datav<(1+datamax-datamin); ++datav)
							{
								for (int modelv=0; modelv<modelGraylevels; ++modelv)
								{
									double pval = Ptablep.getPixelValue(datav,modelv) + pvaloffset;
									fchecksum += pval;
									senddouble( os[i], pval);
								}
							}
						
							os[i].flush();
							pwos[i].println("EndBinary(" + fchecksum + ")");
							pwos[i].flush();
						}
						
						int framestop   = 0;
						int framestart  = 0;
						int last_job = 0;
						for (int i=firstthread; i<use_threads; ++i) //send out Imagedata i=0  may also be send to jobfile
						{
							int threadsum = 0;
							int job = threadVec.get(i).job; //should always be the same job
							int rank = threadVec.get(i).rank;
							int worker = threadVec.get(i).worker;
							if(job != last_job) //reset for each job
							{
								framestop = 0;
								framestart = 0;
								last_job = job;
							}
							
							Idata = jobVec.get(job).Idata;
							IdataSt = Idata.getStack();
							int[] perm = jobVec.get(job).perm;
							int[] sfmaxval = jobVec.get(job).sfmaxval;
							int chsum = 0;
							
							int framecount = jobVec.get(job).portions[rank];
							if(Modus.equals(Modii[1]))
							{
								IJ.log("assigning " + framecount + " frames to worker" + worker);
							}
							
							framestop = framestart + framecount - 1;
							if(skip_stats)
							{
								
								pwos[i].println("updateIdata()");
								pwos[i].flush();
								if(i==0 && debug_level > 3 )
								{
									pwos[use_threads].println("updateIdata()");
									pwos[use_threads].flush();
								}
								int datacount = 3 * (framestop-framestart+1) * modelsize * modelsize / 2;
								pwos[i].println("BeginBinary(" + datacount + ")");
								pwos[i].flush();
								
								if ( 0 == i && debug_level > 3)
								{
									pwos[use_threads].println("BeginBinary(" + datacount + ")");
									pwos[use_threads].flush();
								}
								if(framestop >= perm.length || framestop >= sfmaxval.length)
								{	
									int linenumber = 1 + Thread.currentThread().getStackTrace()[1].getLineNumber(); 
									IJ.log("ERORR: " + getClass().getSimpleName() + " Line" + linenumber + " : too many frames for job " + job);
									framestop = perm.length - 1;
								}
								if (Idata.getType() == ImagePlus.GRAY16)
								{
									for (int frame=framestart; frame<=framestop; ++frame)
									{
										final short[] pixels = (short[])IdataSt.getPixels( perm[frame]);
										for(int ind=0; ind < pixels.length; ++ind)
										{
											if( hex_mask[ind] ) 
											{
												int val = (int)(pixels[ind] & 0xffff);
												send16bitint(os[i],(short)val);
												chsum +=  val;
												if(i==0 && debug_level > 3 )
												{
													send16bitint(os[use_threads],(short)val);
												}
											}		
										}
									}
								}
								else if (Idata.getType() == ImagePlus.GRAY8)
								{
									for (int frame=framestart; frame<=framestop; ++frame)
									{
										final byte[] pixels = (byte[])IdataSt.getPixels( perm[frame]);
										for(int ind=0; ind < pixels.length; ++ind)
										{
											if( hex_mask[ind] ) 
											{
												int val = (int)(pixels[ind] & 0xffff);
												send16bitint(os[i],(short)val);
												chsum += val;
												if(i==0 && debug_level > 3 )
												{
													send16bitint(os[use_threads],(short)val);
												}
											}		
										}
									}
								}
								else if (Idata.getType() == ImagePlus.GRAY32)
								{
									for (int frame=framestart; frame<=framestop; ++frame)
									{
										final float[] pixels = (float[])IdataSt.getPixels( perm[frame]);
										for(int ind=0; ind < pixels.length; ++ind)
										{
											if( hex_mask[ind] ) 
											{
												int val = (int)(pixels[ind]+0.5);
												send16bitint(os[i],(short)val);
												chsum += val;
												if(i==0 && debug_level > 3 )
												{
													send16bitint(os[use_threads],(short)val);
												}
											}		
										}
									}
								}
								else
								{
									IJ.log("Format of Data is not supported");
									saveOptions(options,true);
									return;
								}
								
							}
							else
							{	
								long bytecount = 0;
								long zipcount = 2 * jobVec.get(job).zipImgLen[rank];
								pwos[i].println("updateZipIdata()");
								pwos[i].flush();

								if(0==i && debug_level > 3)
								{
									pwos[use_threads].println("updateZipIdata()");
									pwos[use_threads].flush();
								}

								pwos[i].println("BeginBinary(" + zipcount + ")");
								pwos[i].flush();
								if ( 0 == i && debug_level > 3)
								{
									pwos[use_threads].println("BeginBinary(" +zipcount + ")");
									pwos[use_threads].flush();
								}
								if(framestop >= perm.length || framestop >= sfmaxval.length)
								{	
									IJ.log("ERORR: " + getClass().getSimpleName() + " Line 2168: too many frames for job " + job);
									framestop = perm.length - 1;
								}
								
								//IJ.log("sending "+ framecount  +" frames to worker " + i);//These are correct
								//IJ.log("framestart: " + framestart + " framestop: " + framestop); 
									
								for (int frame=framestart; frame<=framestop; ++frame)
								{
									int sfmax = sfmaxval[ perm[frame]-1 ];
									if (Idata.getType() == ImagePlus.GRAY32)
									{
										//warning: first slice is #1, not #0
										float[] pixels = (float[])IdataSt.getPixels( perm[frame] );
										for(short val = 1; val <= sfmax; ++val)
										{
											short count = 0; //number of occurences of val
											for (short j = 0; j < pixels.length; j++)
											{
												if( (short)(pixels[j] + 0.5) == val ) {++count;}
											}
											if(count > 0)
											{
												send16bitint(os[i],(short)val);
												send16bitint(os[i],(short)count);
												chsum += ( ((int)(val & 0xffff)) + ( (int)(count & 0xffff)) );
												bytecount += 4;
												if(i==0 && debug_level > 3 )
												{
													send16bitint(os[use_threads],(short)val);
													send16bitint(os[use_threads],(short)count);
												}
												for (short j = 0; j < pixels.length; j++)
												{
													if( (short)(pixels[j] + 0.5) == val ) 
													{
														send16bitint(os[i],(short)j);
														if(i==0 && debug_level > 3 )
														{
															send16bitint(os[use_threads],(short)j);
														}	
														chsum += ( (int)(j & 0xffff) );
														bytecount += 2;
													}
												}
											}
										}
										send16bitint(os[i],(short)0);	
										send16bitint(os[i],(short)0);
										bytecount += 4;
										if(i==0 && debug_level > 3 )
										{
											send16bitint(os[use_threads],(short)0);	
											send16bitint(os[use_threads],(short)0);	
										}
									}
									else if (Idata.getType() == ImagePlus.GRAY16)
									{
										short[] pixels = (short[])IdataSt.getPixels( perm[frame]);
										for(short val = 1; val <= sfmax; ++val)
										{
											short count = 0; //number of occurences of val
											for (short j = 0; j < pixels.length; j++)
											{
												if( (short)pixels[j] == val ) {++count;}
											}
											if(count > 0)
											{
												send16bitint(os[i],(short)val);
												send16bitint(os[i],(short)count);
												bytecount += 4;
												chsum += ( ((int)(val & 0xffff)) + ( (int)(count & 0xffff)) );
												if(i==0 && debug_level > 3 )
												{
													send16bitint(os[use_threads],(short)val);
													send16bitint(os[use_threads],(short)count);
												}
												for (short j = 0; j < pixels.length; j++)
												{
													if( (short)pixels[j] == val ) 
													{
														send16bitint(os[i],(short)j);
														if(i==0 && debug_level > 3 )
														{
															send16bitint(os[use_threads],(short)j);
														}	
														bytecount += 2;
														chsum += ( (int)(j & 0xffff) );
													}
												}
											}
										}
										send16bitint(os[i],(short)0);	
										send16bitint(os[i],(short)0);
										bytecount += 4;
										if(i==0 && debug_level > 3 )
										{
											send16bitint(os[use_threads],(short)0);	
											send16bitint(os[use_threads],(short)0);	
										}
									}
									else if (Idata.getType() == ImagePlus.GRAY8)
									{
										byte[] pixels = (byte[])IdataSt.getPixels( perm[frame]);
										for(short val = 1; val <= sfmax; ++val)
										{
											short count = 0; //number of occurences of val
											for (short j = 0; j < pixels.length; j++)
											{
												if( (short)pixels[j] == val ) {++count;}
											}
											if(count > 0)
											{
												send16bitint(os[i],(short)val);
												send16bitint(os[i],(short)count);
												bytecount += 4;
												chsum += ( ((int)(val & 0xffff)) + ( (int)(count & 0xffff)) );
												if(i==0 && debug_level > 3 )
												{
													send16bitint(os[use_threads],(short)val);
													send16bitint(os[use_threads],(short)count);
												}
												for (short j = 0; j < pixels.length; j++)
												{
													if( (short)pixels[j] == val ) 
													{
														send16bitint(os[i],(short)j);
														if(i==0 && debug_level > 3 )
														{
															send16bitint(os[use_threads],(short)j);
														}
														bytecount += 2;
														chsum += ( (int)(j & 0xffff) );
													}
												}
											}
										}
										send16bitint(os[i],(short)0);	
										send16bitint(os[i],(short)0);
										bytecount += 4;
										if(i==0 && debug_level > 3 )
										{
											send16bitint(os[use_threads],(short)0);	
											send16bitint(os[use_threads],(short)0);	
										}
									}
								}
								if(i == use_threads-1)
								{
									if(framestop != (perm.length-1) || framestop != (sfmaxval.length-1))
									{	
										IJ.log("ERORR: " + getClass().getSimpleName() + " Line 2356: too few frames for job " + job);
										return;
									}
								}
								/*
								if(bytecount > zipcount)
								{
									IJ.error("zipped Idata size (" + bytecount + ") exceeds precalculated size ("+ zipcount +")");
									return;
								}
								*/ 
								/*else if(bytecount < zipcount)
								{
									IJ.log("padding zipdata for worker" + i + " with " + ((zipcount-bytecount)/2) + " 0x00 bytes");
								}*/
								
								for ( long padding = bytecount; padding < zipcount; padding += 2 )
								{
									send16bitint(os[i],(short)0);
									if(i==0 && debug_level > 3 )
									{
										send16bitint(os[use_threads],(short)0);			
									}
								}
							}
							
							os[i].flush();
							pwos[i].println("EndBinary(" + chsum + ")");
							pwos[i].flush();
							if ( 0 == i && debug_level > 3)
							{
								os[use_threads].flush();
								pwos[use_threads].println("EndBinary(" + chsum + ")");
								pwos[use_threads].flush();
							}
							//set beginning for next worker:
							framestart = framestop + 1;
							//IJ.log("sent image data to worker" + i);
						}
						
						if(!keep_Idata)
						{	jobVec.get(job_id).Idata.close();}
						jobVec.get(job_id).Ptable.close();
						
						/* *******************************************
						right AFTER receiving Idata a worker would switch from
						jobfile to stdin. Worker will be sleeping now until all jobs are submitted  
						******************************************** */
					}
					
					if(another_job)
					{
						if(modelBrightness > 0)
						{
							IJ.log("WARNING: explicit model scaling and multiple jobs detected");
						}
						
						++job_id; 
						if( batch_run )
						{
							if( (batch_id + 1) < batchpaths.size() ) 
							{	++batch_id;}
							else
							{	
								--job_id;
								System.out.println("ignorring another_job for job " + (batch_id + 1));
								another_job = false; //will quit do {...} while(another_job)
							}
						}
						rank_id = 0;
					}
					if(another_job) //recheck, it may have been ignorred
					{
						IJ.log("appending job " + (job_id + 1) );
						if( Modus.equals(Modii[2]) )
						{
							check_free_nodes = true;
							max_threads = use_threads + loadNodes(path0, false);
						}
					
					}
				}
				while(another_job);
				
				if( job_id > 0 ) //another_job = true
				{
					options = new File(batchpaths.get(0)); //set the master as default
					job_id = 0;
					//retrive original value of job0
					Datatitle = m_Datatitle;
					initthreadnum = m_initthreadnum;
					init_use_frames = m_init_use_frames;
					previous_seed = m_previous_seed;
					seed = m_seed;
					another_job = true;
					sourcepath = m_sourcepath;
					for(int i = 0; i< use_threads; ++i)
					{
						int job =  threadVec.get(i).job;
						int rank = threadVec.get(i).rank;
						int worker = threadVec.get(i).worker;
						IJ.log("Worker" + worker + " has rank" + rank + " on Job" + job);
					}
				}
				
				if( Modus.equals(Modii[1]) )
				{
					saveJobrun(jobrun);
					IJ.log("saved " + jobrun.getPath());
					writeJobStats(jobStats);
					IJ.log("saved " + jobStats.getPath());
					Modus = Modii[2];
					Datatitle = "<none>";
					PtableTitle = "<none>";
					check_free_nodes = false; // jobs are already assigned to nodes!
					saveOptions(options, false);
					IJ.log("saved " + options.getPath());
					Models.close();
					Histimg.close();
					Diffimg.close();
					
					long now = System.currentTimeMillis();
					date.setTime(now);
					IJ.log( dateFormat.format(date) + " " + use_threads + " Job files written in " + format.format((0.001*(now-timeoffset))) + "s" );
					
					saveLog(master_path);
					if(WindowManager.getWindow("Log") != null)
					{
						IJ.selectWindow("Log");
						IJ.run("Close");
					}
					for (int i=0; i<streamnum; ++i)
					{
						pwos[i].close();
						os[i].close();
					}
					sendjobwriter.close();
					sendjobstream.close();
					usenodes_writer.close();
					Process ps1 = rt.exec( "chmod 777 " + sendjob.getPath() );
					ps1.waitFor();
					ps1.destroy();
					continue; //with the batchrun
				}
				//IJ.log("now sending weights and models for " + streamnum + " workers");

				for (int i=0; i<streamnum; ++i) //send out the weights
				{
					try
					{
						pwos[i].println("updateWeights()");
						pwos[i].flush();
						sendWeights( pwos[i], os[i], modelweight);
					}
					catch (Exception e)
					{
						System.out.println("Error: " + threadVec.get(i).command);
						e.printStackTrace();
						return;
					}
					
				}
				//System.out.println("weights OK");
				///TODO implement that in hexworker
				/*
				for (int i=0; i<streamnum; ++i) //send out the inverse relative symmetries
				{
					pwos[i].println("updateSymmetries()");
					pwos[i].flush();
					sendWeights( pwos[i], os[i], matcher.inv_avg_matches);
				}
				*/

				for (int i=0; i<streamnum; ++i) //send out the models
				{
					pwos[i].println("updateModels()");
					pwos[i].flush();
					sendModels(pwos[i], os[i], Models);
				}
				//System.out.println("Models OK");
				for (int i=0; i<streamnum; ++i)
				{
					pwos[i].println("Run(1)");
					pwos[i].flush();
				}
				//System.out.println("Setup OK");
				//read worker data loading report
				//Worker0: Input validated
				for (int i=0; i<use_threads; ++i)
				{
					//print forward messages until end of stream or the expected status line
					do 
					{
						s = reader[i].readLine();
						//if(!s.startsWith("Worker")) System.out.println(s);
						if((s!=null) && (!s.startsWith("Worker"))) 
						{
							System.out.println(threadVec.get(i).command + " " +
							"job: " + threadVec.get(i).job +
							" worker: " + i + " : " + s);
						}
					} 
					while( (s!=null) && (!s.startsWith("Worker") && reader[i].ready() ) );
					if( s == null)
					{
						s = "worker: " + i + " has died, pipe is broken";
						System.out.println(s);	
					}
					IJ.log(s);
				}
				
				{
					long now = System.currentTimeMillis();
					date.setTime(now);
					IJ.log( dateFormat.format(date) + " Data sent in " + format.format((0.001*(now-timeoffset))) + "s" );
					timeoffset = now;
				}
				
				//workers are currently initializing use the time to finalize the nuke file
				nukewriter.close();
				nukestream.close();
				{
					Process ps3 = rt.exec( "chmod 777 " + nuke.getPath() );
					ps3.waitFor();
					ps3.destroy();
				}
				
				double newsum = 0.0;
				double oldsum = 0.0;
				double[] newsums = new double[use_threads];
				double[] oldsums = new double[use_threads];
			
				can_do_stats = (debug_level > 2);
	
				//read initial values from all threads
				for (int i=0; i<use_threads; ++i)
				{
					int worker = threadVec.get(i).worker;
					double seconds = Double.NaN;
					double val = Double.NaN;
					s = reader[i].readLine();
					if(s == null)
					{
						date.setTime(System.currentTimeMillis());
						System.out.println( dateFormat.format(date) );
						IJ.log("Worker " + worker + " has died during initialization aka: " +
								threadVec.get(i).command);
						throw new IllegalStateException("Dead Worker: " + threadVec.get(i).command);
					}
					Scanner sc = new Scanner(s);
					if(sc.hasNextDouble())
					{
						val = sc.nextDouble();
						oldsums[i] = val;
						oldsum += val;
						
						if(sc.hasNextDouble()) seconds = sc.nextDouble();
						if(!Modus.equals(Modii[2]))
						{
							int job =  threadVec.get(i).job;
							int rank = threadVec.get(i).rank;
							jobVec.get(job).performances[rank] = jobVec.get(job).portions[rank]/seconds;
						}
					}
					else
					{
						date.setTime(System.currentTimeMillis());
						System.out.println( dateFormat.format(date) );
						IJ.log("Unexpected reply from Worker:" + worker + " aka " + threadVec.get(i).command);
						IJ.log(s);
						throw new IllegalStateException("Unexpected reply: " + s + " from" + threadVec.get(i).command);
					}
					System.out.println( threadVec.get(i).command + "\tworker" + worker + ": " + format.format(val) + " in " + format.format(seconds) + "s");
				}
				
				for(int m=0; m < modelnum; ++m)
				{	cases[m] = 0;}				
				for (int i = 0; i < use_threads; ++i) //Dont send initial Pollcases to job.bin
				{
					pwos[i].println("PollCases()");
					pwos[i].flush();
				}	
				for (int i = 0; i < use_threads; ++i) //collect initial cases
				{	
					for(int m=0; (m < modelnum) ; ++m)
					{
						cases[m] += Integer.parseInt( reader[i].readLine() );
					}
				}
				if(!use_weights)
				{
					//System.out.println("cases -> weights");
					for(int m=0; (m < modelnum) ; ++m)
					{
						modelweightupd[m] = (double)cases[m];
					}
					normalize(modelweightupd);
					for(int m=0; (m < modelnum) ; ++m)
					{
						modelweight[m] = modelweightupd[m];
					}
					saveWeights(master_path);
				}
				/*
				if(!Modus.equals(Modii[2]))
				{
					writePerformance(master_path);
					System.out.println("saved Performance.txt");	
				}
				*/
				IJ.log("Pval = " + format.format(oldsum));
				//System.out.println("Pval =\t" + format.format(oldsum));
				{
					long now = System.currentTimeMillis();
					date.setTime(now);
					IJ.log( dateFormat.format(date) + " Initialization done in " + format.format((0.001*(now-timeoffset))) + "s");
					timeoffset = now;
				}
				Pval = oldsum;
				highestPval = Pval;
				//////////////////////////////////////////
				//here comes the main loop
				/////////////////////////////

				//debugmsg("Elapsed time:",(System.currentTimeMillis()-timeoffset)/1000.);

				//the list of changes
				int chlistmax = 2*(blocksize+2)*(blocksize+2)*blockdepth; //safely includes blobbed zone
				int chlistlen = 1;
				int[] chx  = new int[chlistmax];
				int[] chy  = new int[chlistmax];
				int[] chm  = new int[chlistmax];
				int[] newv = new int[chlistmax];
				int[] oldv = new int[chlistmax];
				//double avg_guess = 0.0;

				int loopcounter = 0;
				int loopmax = report_interval;
				request_cache_report = (debug_level > 0);

				//Adjust the same pixel in parallel for all models

				incrs = new int[blockvolume];
				testval = new int[blockvolume];    //the current guess
				oldval = new int[blockvolume];  //the last guess
				oldvalinit = new int[blockvolume]; //the starting point
				bestval = new int[blockvolume]; //the best guess so far
				secondbestval = new int[blockvolume]; //the secondbest guess so far
				thirdbestval = new int[blockvolume]; //the secondbest guess so far
				avg_blockval = new int[blockvolume]; 
				bestPval = new double[blockvolume]; //the highest Pval
				secondbestPval = new double[blockvolume]; //the secondhighest Pval
				thirdbestPval = new double[blockvolume]; //the thirdhighest Pval
				
				double next_factor = 1.0;
				double last_factor = 1.0;
				
				//default are complete models
				
				
				if(spiral)
				{
					spiralcountmax = 1;
					spiralDelta = 1;
					spiralcount = -1;
					spiraldirection = 0;
					homex = modelsize/2-blocksize;
					homey = modelsize/2-blocksize;
					if(flip_model) {homey -= blocksize;}
					else {homex -= blocksize;}
				}
				
				if(wrapnum == 0 && (!do_full_EM))
				{ --wrapnum;}
				
				
				homem = 0;
				homex = 0;
				homey = 0;
				
				
				if(shuffle_update)
				{
					int blockind = get_and_shuffle(++blockid, permutation[homem]);
					//permutator.getIndex(++blockid);
					int dim = modelsize/blocksize;
					homex = blocksize * (blockind % dim) + ( ((wrapnum % 2) == 1)?blocksize/2:0 );
					homey = blocksize * (blockind / dim) + ( ((wrapnum % 2) == 1)?blocksize/2:0 );
				}
				blockx = 0;
				blocky = 0;
				blockm = 0;
				
				strategy = Strategy.PIXELADJ; //that should be the only strategy for the first complete wrap
				initadjustmentBlock(Models, Activity);
				
				char[] pmsymb = {'-', '0', '+', '!'};
				weightupdnum = 0 ;
				
				int progstep = loopmax / 100;
				if(progstep == 0)
				{	progstep = 1;}
				
				boolean clear_updates = true;
				boolean wraps_completed = false;
				boolean get_new_weights = false;
				boolean recent_em = false;
				mutex = 1;
				
				do //mainloop
				{
					
					boolean tentative_acceptance = false;
					if(recent_em)
					{
						reportWrap();
						recent_em = false;
						clear_updates = true;
						if( !wraps_completed && ( target_wraps > -1) && ( wrapnum >= target_wraps) )
						{	
							IJ.log("reached target for number of wraps by EM step, polishing weights ... ");
							wraps_completed = true;
							weight_trials = 100*modelnum;
						}
					}
					
					boolean accept_always = false;
					do //strategyloop
					{
						chlistlen = 0;
						cmsok = 0; //if we have a new request
						//choose a new strategy and/or block	
							
						if( (renew==1) && (ongoing_Adj == -1) && (weight_trials < 1) && (move_trials <1) ) //no more active pixels in the current block
						{
								strategy = Strategy.PIXELADJ;
								boolean advBlock = (advanceBlock() == 1);
								clear_updates |= advBlock; //may trigger reportWrap() which resets weight_trials
								if( !wraps_completed && ( target_wraps > -1) && ( wrapnum >= target_wraps) )
								{	
									IJ.log("reached target for number of wraps");
									wraps_completed = true;
									if(use_weights && (!const_weights) && (!aggressive_weights))
									{	
										weight_trials = 100*modelnum;
										IJ.log("polishing weights: " + weight_trials);
									}
								}
								
								if( advBlock && ( (new_blocksize != blocksize) || (new_blockdepth != blockdepth) ) )
								{
									//nasty code duplicate from above
									blocksize = new_blocksize;
									initblocksize = blocksize;
									blockdepth = new_blockdepth;
									initblockdepth = blockdepth;
									blockarea = blocksize * blocksize;
									blockvolume = blockarea * blockdepth;	
									blockid = -1;
									if(shuffle_update)
									{   //Third copy of that code, NEEDS REFCATORING see advanceBlock()
										int num_blocks = (modelsize * modelsize) / blockarea;
										permutation = new int[modelnum][num_blocks];
										init_perm(permutation);
										//permutator = new Permutator(num_blocks,num_blocks/2); //first block is in center
										int blockind = get_and_shuffle( ++blockid, permutation[homem]); //actually always 0
										//permutator.getIndex(++blockid);
										int dim = modelsize/blocksize;
										homex = blocksize * (blockind % dim) + ( ((wrapnum % 2) == 1)?blocksize/2:0 );
										homey = blocksize * (blockind / dim) + ( ((wrapnum % 2) == 1)?blocksize/2:0 );	
									}
									chlistmax = 2*(blocksize+2)*(blocksize+2)*blockdepth;
									chx  = new int[chlistmax];
									chy  = new int[chlistmax];
									chm  = new int[chlistmax];
									newv = new int[chlistmax];
									oldv = new int[chlistmax];
									
									incrs = new int[blockvolume];
									testval = new int[blockvolume]; //the current guess
									oldval = new int[blockvolume];  //the last guess
									oldvalinit = new int[blockvolume]; //the starting point
									bestval = new int[blockvolume]; //the best guess so far
									secondbestval = new int[blockvolume]; //the secondbest guess so far
									thirdbestval = new int[blockvolume]; //the secondbest guess so far
									avg_blockval = new int[blockvolume];
									bestPval = new double[blockvolume]; //the highest Pval
									secondbestPval = new double[blockvolume]; //the secondhighest Pval
									thirdbestPval = new double[blockvolume]; //the thirdhighest Pval
								}
								initadjustmentBlock(Models, Activity);
						}
						//moves shall be completed before weights start
						//at the moment we are nowhere explicitely requiering moves
						//the can however be issue by maketestadjustment
						if( (move_trials > 0) ) 
						{
							--move_trials;
							strategy = Strategy.PIX_MOVES;
						}
						//weights have to be completed before pixel adjustment can proceed
						//weight_trials could be issues by AdvanceBlock
						else if( weight_trials > 0)
						{
							--weight_trials;
							strategy = Strategy.WEIGHTADJ;
							//TODO check if we could spare this call if tell_use_weights is true
							if(use_weights && init_weights/*!tell_use_weights*/)
							{	
								initWeightAdjustment(true);
								init_weights = false;
							}						
						}
						else //resume pixeladjustments or EM
						{
							
							if(strategy != Strategy.PIXELADJ)
							{
								if( wraps_completed || !do_full_EM  /*|| (wrapnum<=modelnum/multiActivations)*/ )
								{
									if( collect_Histos && (strategy==Strategy.WEIGHTADJ) ) //&& (!wraps_completed)
									{
										strategy = Strategy.POLL_HISTO;
									}
									else
									{
										strategy = Strategy.PIXELADJ;
									}
								}
								else
								{
									strategy = Strategy.TOTAL_EM;
									recent_em = true;
								}		
							}
						}	
						
						
						int i3d = blockx + blocky * blocksize + blockm * blockarea;
						switch(strategy) //follow the strategy plans
						{
							case WEIGHTADJ:
								renew = 0; //only once at a time
								cmsok = 1;
								if(use_weights && !tell_use_weights)
								{	weightAjustment();}
								break;
							case PIXELADJ:
								renew = 0;
								//stay on Pixel until it is converged
								if((incrs[i3d]!=-100) && (correlated_optimization || (incrs[i3d] == 3) )) //|| force_wrap
								{
									if((advancePixel()==1))
									{
										
										if (advanceModel() == 0) //next model
										{
											initadjustmentBlock(Models, Activity);
											
													
										}
										else //run over all models finished
										{
											renew = 1;
											cmsok = 0;
											ongoing_Adj = -1;
											break;
										}
									}
								}
								if(ongoing_Adj > 0)
								{
									short[] pixels= (short[])ModelsSt.getPixels(((homem + blockm + modelnum) % modelnum) + 1);
									int bx = (homex + blockx + modelsize) % modelsize;
									int by = (homey + blocky + modelsize) % modelsize;
									int currentval = pixels[ bx + by * modelsize ];
									boolean commit_EM = false;
									if(incrs[i3d] == -1 && !wraps_completed)
									{
										if(do_pix_EM && force_pix_EM)
										{
											commit_EM = true; //always do pix_EM if it is forced i.e. always accepted
										}
										else if (do_pix_EM)
										{
											if(random.nextDouble() <= 0.01)
											{	commit_EM = true;}
										}
										
									}
									
									if( !commit_EM )
									{
										//might switch strategy to PIX_MOVES
										chlistlen = mktestadjustment3(chx,chy,chm,newv,currentval,chlistmax,chlistlen);
										//tentative_acceptance = random.nextBoolean();  //what to do if workers report NCH
										if(blobs)
										{
											chlistlen = blobify_changes(Models, chlistlen, chx, chy, chm, newv, oldv);
										}
										int mbind = blockx + blocky * blocksize + blockm * blockarea;
										tentative_acceptance = ( Math.abs(newv[0] - avg_blockval[mbind]) <=
																 Math.abs(oldv[0]-avg_blockval[mbind])	);
									}
									else
									{
										strategy = Strategy.PIXEL_EM;
										incrs[i3d] = -100;
										//System.out.println("choosing strategy = PIXEL_EM");
									}
								}
								break;
								case POLL_HISTO:
									renew = 0;
									cmsok = 1;
								break;
								case PIX_MOVES: //explicit case
									renew = 1; //only once at a time
									cmsok = 1;
									chlistlen = switch_pixels(Models, chx, chy, chm, newv, oldv);
									++dbg_moves;
									if(chlistlen==0)
									{
										System.out.println("Warning: switch_pixels did not succeed");
										cmsok = 0; //pick another strategy
									}
									if(blobs)
									{
										chlistlen = blobify_changes(Models, chlistlen, chx, chy, chm, newv, oldv);
									}
								break;
						}		
						
						switch(strategy) //may have been changed depdening on incrs;
						{	
							
							case PIXEL_EM:
								chx[0] = (homex + blockx + modelsize) % modelsize;
								chy[0] = (homey + blocky + modelsize) % modelsize;
								//falltrough
							case MODEL_EM:
								chm[0] = (homem + blockm + modelnum) % modelnum;
								
								//falltrough	
							case TOTAL_EM:
								chlistlen = 1;
								//fallthrough 
							case PIX_MOVES:
								//these strategies will be commit once and then reevaluated
								renew = 1;
								cmsok = 1;
								
								
								break;
							case POLL_HISTO:
								chlistlen = 0;
								break;
							
						}

						



						if( ( strategy == Strategy.PIXELADJ ) || (strategy == Strategy.PIX_MOVES) )//check for valid requests
						{
							for (int q=0; q<chlistlen; ++q) //fill in oldv[q] and check != newv[q]
							{
								short[] pixels= (short[])ModelsSt.getPixels(chm[q]+1);
								oldv[q] = pixels[ chx[q] + modelsize*chy[q] ];
								int ind3d = blockx + blocky * blocksize + blockm * blockarea;
								if(oldv[q] == newv[q])
								{
									System.out.printf("WARNING, useless update request reply:%d duplicate:%d\n", reply, oldv[q]);
									System.out.print("Strategy " + getStrategyTag());
									System.out.print(" incrs " + incrs[ind3d]);
									System.out.println(" home x,y,m " + homex + "," + homey + "," + homem + " block x,y,m " + blockx + "," + blocky + "," + blockm);
									System.out.print("olvalinit: " + oldvalinit[ind3d]);
									System.out.print(" oldval: " + oldval[ind3d]);
									System.out.print(" bestval: " + bestval[ind3d]);
									System.out.print(" secondbestval: " + secondbestval[ind3d]);
									System.out.print(" thirdbestval: " + thirdbestval[ind3d]);
									System.out.println(" avg_blockval: " + avg_blockval[ind3d]);
									cmsok = 0; //cancel the plan
								}
								if( (modelnum > 1) && (wrapnum < 1) && (chm[q] == 0) && (modelnum>multiActivations) )
								{
									System.out.println("WARNING, premature update for model0 x,y " + chx[q] + "," + chy[q]);
									incrs[ind3d] = 3;//declare the pixel as converged for this wrap
									--ongoing_Adj;
									cmsok = 0; //cancel the plan
								}
								if(chx[q] < 0 || chx[q] > modelsize-1 ||
									chy[q] < 0 || chy[q] > modelsize-1
									|| newv[q] < 0 || newv[q] > modelGraylevels - 1
									|| chm[q] < 0 || chm[q] > modelnum - 1)
								{
									date.setTime(System.currentTimeMillis());
									System.out.println( dateFormat.format(date) );
									System.out.print("ERROR, invalid changes detected ");
									System.out.println("chx " + chx[q] + ", chy " + chy[q] + ", chm " + chm[q] + ", newval " + newv[q] + ", oldval: " + oldv[q]);
									System.out.print("Strategy " + getStrategyTag());
									System.out.println("chlistlen: " + chlistlen);
									System.out.print(" incrs " + incrs[ind3d]);
									System.out.print(" homem,x,y " + homem + "," + homex + "," + homey);
									System.out.println(" blockm,x,y " + blockm + "," + blockx + "," + blocky);
									cmsok = 0; //cancel the plan
									throw new IllegalStateException("coordinates for changes are out of range.");
								}
							} //DONT use oldv for applying changes
						}
						
					
						//we are at the last weight step AND it has converged
						boolean commit_autosave = ( request_autosave &&
													(weight_trials < 1) &&
													strategy != Strategy.WEIGHTADJ );						 
						if(commit_autosave)
						{	//TODO turn this duplicate code a function
							request_autosave = false; //is set true by wrapreport		
							System.out.println(master_path);
							saveOptions(options, false);
							System.out.println("saved " + options.getName());
							saveWeights(master_path);
							
							if( can_do_stats )
							{	
								writeSignificance(master_path, Models);
								System.out.println("saved " + "Significance.txt");
							}
							recent_save = true;
							saveAllImgs( master_path );
							try
							{
								Process cpPval = rt.exec("cp -f " + master_path + "Pval.txt " + master_path + "Pval.bak"); 
								cpPval.waitFor();
								cpPval.destroy();
								System.out.println("Pval.txt -> Pval.bak");
							}
							catch( Exception e)
							{
								e.printStackTrace();
								//actually no serious problem
							}
							
							recent_img_save = true;
						}
						
						if (cmsok==0) { renew = 1;}
					
					}
					while (cmsok==0);
					
					int ind3d = blockx + blocky * blocksize + blockm * blockarea;
					if(request_cache_report) //is set true inside reportWrap()
					{
						request_cache_report = false;
						for (int i = 0; i < streamnum; ++i) 
						{
							if(can_write[i])
							{
								pwos[i].println("CacheReport(" + wrapnum + ")");
								pwos[i].flush();
							}
						}
					}
					//boolean force_no_changes = false; 
					if(strategy == Strategy.PIXELADJ)
					{
						if(tell_use_max)
						{
							//force_no_changes = true;
							accept_always = true;
							tell_use_max = false;
							if(shuffle_update)
							{	--blockid;}
							chlistlen = 0;
							IJ.log("use max: " + use_max);
							System.out.println("use max: " + use_max);
							for (int i = 0; i < streamnum; ++i) 
							{
								pwos[i].println("use_max(" + (use_max?1:0)  + ")");
								pwos[i].flush();
							}	
						}
						if(tell_use_weights) //might happen if reading jobfiles with wrong settings
						{
							tell_use_weights = false;
							IJ.log("use weights: " + use_weights + ( (use_weights && const_weights)?" const":""));
							System.out.println("use weights: " + use_weights + ( (use_weights && const_weights)?" const":""));
							for (int i = 0; i < streamnum; ++i) 
							{
								pwos[i].println("use_weights(" + (use_weights?1:0)  + ")");
							}
							//Polling should not be necessary since jobrun requires Weights.txt
						}
						
						next_factor = make_nextModel(Models, chlistlen, chx, chy, chm, newv, oldv);
						last_factor = sym_factors[chm[0]]; //backup copy
						sym_factors[chm[0]] = next_factor;
					}
					
					if( (clear_updates || wraps_completed) )
					{
						//System.out.println("polling cases");
						for(int m=0; m < modelnum; ++m)
						{	cases[m] = 0;}
							
						for (int i = 0; i < streamnum; ++i)
						{
							pwos[i].println("PollCases()");
							pwos[i].flush();
						}	
						for(int i = 0; i < use_threads; ++i)
						{
							for(int m=0; (m < modelnum) ; ++m)
							{
								cases[m] += Integer.parseInt( reader[i].readLine() );
							}
						}	
					} 
					
					
					if( (strategy == Strategy.WEIGHTADJ) )
					{
						if(tell_use_weights)
						{
							IJ.log("use weights: " + use_weights);
							System.out.println("use weights: " + use_weights + ( (use_weights && const_weights)?" const":"") );
							tell_use_weights = false;
							for (int i = 0; i < streamnum; ++i) 
							{
								pwos[i].println("use_weights(" + (use_weights?1:0)  + ")");
							}
							get_new_weights = true; //we want new weights AND they shall be accepted
							
						}
						
						if( weight_trials == 0 ) 
						{	
							
							if(get_new_weights || !use_weights)
							{
								System.out.println("cases -> weights");
								for(int m=0; (m < modelnum) ; ++m)
								{
									modelweightupd[m] = (double)cases[m];
								}
								normalize(modelweightupd);
								initWeightAdjustment(false);
								get_new_weights = false;
								accept_always = true;	
								renew = 1;
							}	
						}
					}
					
					
					
					/*if(strategy == Strategy.TOTAL_EM)
					{	
						System.out.println("Broadcasting request for total EXPECTATION MAXIMATION now!");
					}*/
					
					for (int i = 0; i < streamnum; ++i)
					{
						boolean commit_run = true;
						
						
						switch (strategy) //broadcast the respective commands
						{
							case WEIGHTADJ:
								pwos[i].println("updateWeights()");
								pwos[i].flush();
								sendWeights(pwos[i], os[i], modelweightupd); //scaled tentative weights
								break;
							case PIX_MOVES:
							//fallthrough
							case PIXELADJ:
								if( (commit_run = !wraps_completed) )
								{
									if(sym_scaling)
									{
										pwos[i].println("updateWeights()");
										pwos[i].flush();
										sendWeights(pwos[i], os[i], modelweight); //scaled current weights
									}
									
									for (int q = 0; q < chlistlen; ++q)
									{
										pwos[i].println("ModelPixel("+chm[q]+","+chx[q]+","+chy[q]+","+newv[q]+")"); 
									}
								}
								break;
							case PIXEL_EM:
								if( (commit_run = !wraps_completed) )
								{
									//if(i==0)	System.out.println("Pixel_EM("+chm[0]+","+chx[0]+","+chy[0]+")");
									for (int q = 0; q < chlistlen; ++q)
									{
										pwos[i].println("Pixel_EM("+chm[q]+","+chx[q]+","+chy[q]+")"); 
									}
								}
								break;		
							case MODEL_EM:
								if( (commit_run = !wraps_completed) )
								{
									for (int q = 0; q < chlistlen; ++q)
									{
										//ranges from first to last model are supported
										pwos[i].println("Model_EM("+chm[q]+","+chm[q]+")");
									}
								}
								break;
								
							case TOTAL_EM:
								if( (commit_run = !wraps_completed) )
								{
									for (int q = 0; q < chlistlen; ++q)
									{
										pwos[i].println("Model_EM("+ 0 + "," + (modelnum-1) + ")");
									}
								}
								break;
							case POLL_HISTO:
								commit_run = true;
								pwos[i].println("updateHisto");
								break;	
							default:
								date.setTime(System.currentTimeMillis());
								System.out.println( dateFormat.format(date) );
								System.out.println("default in STRATEGY while braodcasting commands");
								throw new IllegalStateException("undefined STRATEGY while broadcasting");
						}
						
						if(commit_run)
						{
							pwos[i].println("Run(1)");
							pwos[i].flush();
						}
					}
					fresh_symmetries = false;
					//workes are busy now, do the gui stuff now
					pending_keep = !( (strategy == Strategy.PIXELADJ || strategy == Strategy.PIX_MOVES) && wraps_completed);
					switch(strategy) //these never require a "keep" acknowledgement
					{
						case PIXEL_EM:
						case MODEL_EM:
						case TOTAL_EM:
						case POLL_HISTO:
							pending_keep = false;
						break;
					}
					
					Activity.updateAndDraw();
					
					/*if((loopcounter % progstep) == 0)
					{
						IJ.showProgress(loopcounter, loopmax);
					}*/
					if(strategy == Strategy.PIXELADJ || strategy == Strategy.TOTAL_EM ||
					   strategy == Strategy.MODEL_EM || strategy == Strategy.POLL_HISTO ||
					   strategy == Strategy.PIX_MOVES)
					{
						//keep the image during the weight steps
						if(clear_updates || wraps_completed) 
						{
							clear_updates = false;
							for(int m = 1; m <= modelnum; ++m)
							{
								update_pixels = (short[])UpdatesSt.getPixels(m);
								for(int i = 0; i < update_pixels.length; ++i)
								{
									update_pixels[i] = 0;
								}
							}
							int fakenum = (wrapnum > -1)?wrapnum:0;
							System.out.print("Wrap: " + fakenum + "\tWeights:");
							if(use_weights)
							{	System.out.print(" active" + ( const_weights?" (const)":"") + (aggressive_weights?" (aggressive)":"") ); }
							System.out.println();
							
							for (int v = 0; v < modelnum; ++v)
							{
								System.out.print("" + format.format(modelweightupd[v]));
								if (dead[v] == ModelState.FINISHED)
								{	System.out.print("'");}
								//if(v == weightupdnum)
								//{	System.out.print("_");}
								if ( (v + 1) % weight_cols == 0 || (v == modelnum - 1) )
								{	System.out.println();}
								else
								{	System.out.print("\t");}
							}
							if(sym_scaling)
							{
								System.out.println("Wrap: " + fakenum + "\tweight scaling: " + sym_scaling + "\tmedian filtering: " + median_filtering);
								for (int v = 0; v < modelnum; ++v)
								{
									System.out.print("" + format.format(v==chm[0]?last_factor:sym_factors[v]));
									if (dead[v] == ModelState.FINISHED)
									{	System.out.print("'");}
									if ( (v + 1) % weight_cols == 0 || (v == modelnum - 1) )
									{	System.out.println();}
									else
									{	System.out.print("\t");}
								}
							}
							if(wrapnum > -1)
							{
								System.out.print("Wrap: " + fakenum + "\tCases:");
								if(!use_weights)
								{	System.out.print(" active"); }
								System.out.println();
								for (int v = 0; v < modelnum; ++v)
								{
									System.out.print("" + cases[v] );
									if (dead[v] == ModelState.FINISHED)
									{	System.out.print("'");}
									if ( (v + 1) % weight_cols == 0 || (v == modelnum - 1) )
									{	System.out.println();}
									else
									{	System.out.print("\t");}
								}
							}
							if(use_weights || wrapnum == -1 )
							{
								System.out.println("Pval = " + format.format(oldsum));
							}
							if(use_weights && (!const_weights) && (wrapnum > 0) )
							{
								if(dbg_weights > 0)
								{	IJ.log("weights : " + dbg_weights + "    Pval = " + format.format(oldsum));}
								dbg_weights = 0;
							}
						}
						
						//show current change on Updates Image
						if(pending_keep)
						{
							for (int q = 0; q < chlistlen; ++q)
							{
								update_pixels = (short[])UpdatesSt.getPixels(chm[q]+1);
								update_pixels[ chx[q] + modelsize * chy[q] ] += 1;
								/*
								///DEBUG Section somehow the block mover got stuck
								if(update_pixels[ chx[q] + modelsize * chy[q] ] > 100)
								{
									System.out.println("massive updates detected");
									System.out.println("strategy: " + getStrategyTag());
									System.out.println("at (" + chx[q] + "," + chy[q] + ") in model: " + chm[q] );
									System.out.println("renew: " + renew);
									System.out.println("cmsok: " + cmsok);
									System.out.println("ongoing_Adj: " + ongoing_Adj);
									System.out.println("weight_trials: " + weight_trials);
									System.out.println("move_trials: " + move_trials);
									
									
									return;
								}
								*/ 
							}
							Updates.updateAndDraw();
							Models.updateAndDraw();
							Diffimg.updateAndDraw();
						}
					}
					
					if (writehist) //update model history if needed
					{
						//if(wrapnum < 1)
						//{	System.out.println("WARNING, inavlid request for writing to Model History");}
						writehist = false;
						int histlength = Histimg.getStackSize();  //1 is first one
						ImageStack HistimgSt = Histimg.getStack();
						float[] histpix = new float[Histimg.getWidth() * Histimg.getHeight()];
						float[] histpixprev = (float[])HistimgSt.getPixels(histlength);
						int modelarea = modelsize * modelsize;
						for (int m = 0; m < modelnum; ++m)
						{
							short[] pixels = (short[])ModelsSt.getPixels(m+1);
							short[] diffpix = (short[])DiffSt.getPixels(m+1);
							for(int ind = 0; ind < modelarea; ++ind)
							{
								int x = ind % modelsize;
								int y = ind / modelsize;
								float val = /*(float)modelmin0+*/((float)pixels[x+modelsize*y])*
									((float)(modelmax0-modelmin0))/((float)modelGraylevels);
									
								int intdiff = ( diffpix[x+modelsize*y] & 0xffff );
								double diffval = ( (double)(intdiff - 32768) ) * 
												(modelmax0-modelmin0) / 
												( (double)modelGraylevels ); //actual intensity
								histpix[x+m*(modelsize+1)+histwidth*y] = val;
								histpix[x+m*(modelsize+1)+histwidth*(y+modelsize+1)] = (float)diffval;
								histpix[x+m*(modelsize+1)+histwidth*(y+2*modelsize+2)]
									=val-histpixprev[x+m*(modelsize+1)+histwidth*y];
							}
						}
						HistimgSt.addSlice("Wrap" + wrapnum, histpix);					
						Histimg.setStack(HistimgSt);
					}

					if(  wraps_completed && (strategy == Strategy.PIXELADJ) ) //break from mainloop
					{	break; }


					//print the first two chars of console report
					if (loopcounter > (loopmax-report_num + 1) )
					{

						switch(strategy)
						{
							case PIXELADJ:
								//all pixel test are assumed to happen on the same model, the clients validate that anyways
								System.out.print("" + (incrs[ind3d]) + "M" + chm[0] + " " ); 
							break;
							case WEIGHTADJ:
								System.out.print("" + weight_incrs + "W  " );
							break;
							case PIX_MOVES:
								System.out.print("SW"+chm[0] );
							break;
							case PIXEL_EM:
								//System.out.print("Pixel_EM");
								break;
							case MODEL_EM:
								//System.out.print("Model_EM");
								break;
							case TOTAL_EM:
								//System.out.print("Total_EM");
								break;
							case POLL_HISTO:
								//System.out.print("updateHisto");
								break;			
							default:
								date.setTime(System.currentTimeMillis());
								System.out.println( dateFormat.format(date) );
								throw new IllegalStateException("undefined STRATEGY while reporting to console");
						}
					}
					boolean force_update = false;
					// force update may be triggered by user at bottom of the loop
					do //typically a one time loop
					{
						int fm = 0;
						int lm = (modelnum-1);
						switch(strategy)
						{
							case PIX_MOVES:
							//falltrough
							case WEIGHTADJ:
							//falltrough
							case PIXELADJ:
							
								newsum = 0.0;
								for (int i=0; i<use_threads; ++i) //read Pvalues from workers and update newsums
								{
									int worker = threadVec.get(i).worker;
									s = reader[i].readLine();

									if (s == null) // report error and return
									{
										date.setTime(System.currentTimeMillis());
										System.out.println( dateFormat.format(date) );
										System.out.println("NullString worker: " + worker + getStrategyTag() + (incrs[ind3d]) + 
										" aka: " + threadVec.get(i).command);
										IJ.log("NullString worker: " + worker + getStrategyTag() + (incrs[ind3d]) + 
										" aka: " + threadVec.get(i).command);
										throw new IllegalStateException("Unexpected reply");
									}
									double thisAmount = 0.0;
									Scanner sc = new Scanner(s);
									String tag = sc.next();
									if(tag.equals("PVAL") && sc.hasNextDouble() )
									{
										thisAmount = sc.nextDouble();	
									}	
									else if(tag.equals("NCH") )
									{
										thisAmount = oldsums[i];
										tentative_acceptance &= use_beamprofile;
											
									}
									else //report error and return
									{
										date.setTime(System.currentTimeMillis());
										System.out.println( dateFormat.format(date) );
										IJ.log("Unexpected reply (PVAL) from Worker:" + worker + " aka: " + threadVec.get(i).command);
										IJ.log(s);
										System.out.println("received: " + s);
										System.out.println("from: " + threadVec.get(i).command);
										throw new IllegalStateException("Unexpected reply");
									}

									newsum += thisAmount;
									newsums[i] = thisAmount;
									if ( loopcounter > (loopmax-report_num + 1) )
									{
										int compthr = Double.compare(newsums[i], oldsums[i]);
										System.out.print(pmsymb[1+compthr]);
									}
								}

								reply = -1;
								reply = (accept_always) ? 1 : Double.compare(newsum, oldsum);
								
								if(incrs[ind3d]==-200) //we are testing an EM value
								{
									
									if(force_pix_EM)
									{	reply = 1;}
									else if(debug_level > 1)
									{	
										System.out.println( ( (reply==1)?"accepted":"rejected") + " with dPval = " + format.format(newsum-oldsum) );
									}
								}
								
								Pval = newsum; //global copy
								if(force_update)
								{
									reply = 1;
									System.out.println("old Pval= " + format.format(oldsum));
									System.out.println("new Pval= " + format.format(newsum));
									IJ.log("old Pval= " + format.format(oldsum));
									IJ.log("new Pval= " + format.format(newsum));
								}
								
								if( (reply==0) && ( (tentative_acceptance) /*|| (testing_avg && (incrs[mbind]==0) )*/ ) )
								{
									/*
									int mbind = blockx + blocky * blocksize + blockm * blockarea;
									
									if(testing_avg && (incrs[mbind]==0) )
									{
										System.out.println("accepting minor averaging");
									}
									else
									{
										System.out.println("accepting smoothing minor modelchanges");
									}
									*/ 
									reply=1;
								}
								/*else if(reply==0)
								{
									System.out.println("rejecting roughing minor modelchanges");
								}*/
								//(sending integer, 1=yes, -1=no (0:equal))
								for (int i = 0; i < streamnum; ++i)
								{
									pwos[i].println("Keep(" + reply + ")");
									pwos[i].flush();
								}
								pending_keep = false;
								if(!force_update && (reply == 1) && ((strategy == Strategy.PIXELADJ) || (strategy == Strategy.PIX_MOVES) ))
								{
									///HMMM this also increases with blobbed chlistlen, maybe good or bad?
									/// rather bad way to many updates in every run
									//make all pixels around an intensity swap active
									final int cores = (strategy == Strategy.PIX_MOVES)?chlistlen:1;
									//The entire state of the active adjustment block
									//is about to be messed up. just set generous activity and then move on
									if(strategy == Strategy.PIX_MOVES) 
									{
										//System.out.print("$" + chm[0]);
										//System.out.flush();
										++dbg_good_moves;
										ongoing_Adj = -1;
										renew = 1;
										cmsok = 0;
										move_trials = 0;
										weight_trials = aggressive_weights?4*modelnum:0; 
									}
									
									for (int q = 0; q < cores; ++q)
									{
										activity_pixels = (short[])ActivitySt.getPixels(chm[q]+1);
										int val = activity_pixels[ chx[q] + modelsize * chy[q] ] & 0xffff;
										activity_pixels[ chx[q] + modelsize * chy[q] ] = (short)(val + 96); //reward for change
										int x = chx[q] - chy[q]/2;
										int z = chy[q];
										int y = -x -z;
										//actually inefficent looping but totally ok for radius = 1
										if( (wrapnum >= modelnum/multiActivations) || 
											(distToCenter(chx[q], chy[q]) < modelactive/2) ) // dont spill activity to silent ring during launches
										{
											for(int dx = -1; dx <= 1; ++dx)
											{
												for(int dz = -1; dz <= 1; ++dz)
												{
													int dy = -dx - dz;
													if(Math.abs(dx)+Math.abs(dy)+Math.abs(dz) == 2)
													{
														int nx = x + dx;
														int nz = z + dz;
														int ny = -nx -nz;
														int[] cube = {nx-cx,ny-cy,nz-cz};
														hexagonaltrap(cube, modelsize/2);
														nx = cube[0] + cx;
														ny = cube[1] + cy;
														nz = cube[2] + cz;
														int nq = nx + nz/2;
														int nr = nz;
														int ind = nq + nr * modelsize;
														int val2 = activity_pixels[ ind ] & 0xffff;
														activity_pixels[ ind ] = (short)(val2 + 16); //reward for neighboring change
														
													}	
												}
											}
										} 	
									}
								}
								
								if (loopcounter > (loopmax - report_num + 1) && !force_update)
								{
									int run = loopcounter - loopmax + report_num - 1; //starts with 1
									System.out.print("|"+pmsymb[1+reply]+" "+chlistlen);
									if( (run % report_cols == 0) || (run == report_num) )
									{	System.out.println();}
									else
									{	System.out.print("\t");}
								}

								if (loopcounter > loopmax && !force_update) //print Pval
								{
									long now = System.currentTimeMillis();
									System.out.println("Pval= " + format.format(oldsum) + "\tTests/s: " + 
									format.format( (double)(pixel_tests + dbg_weights)/(0.001*(now-timeoffset)) ) );
									IJ.log("Pval= " + format.format(oldsum) + "  Tests/s: " + 
									format.format( (double)(pixel_tests + dbg_weights)/(0.001*(now-timeoffset)) ) );
									loopcounter = 0;
									if(initbenchmark) //ok we quit after first Pval
									{
										date.setTime(now);
										IJ.log( dateFormat.format(date) + " reached second Pval in " + format.format((0.001*(now-timeoffset))) + "s");
										IJ.log("Tests/s: " + format.format( (double)(pixel_tests + dbg_weights)/(0.001*(now-timeoffset)) ) );
										//timeoffset = now;
										domainloop = false;
										break; //from the mainloop
									}
								}
								//fallthrough
							case CONFIRMED_PIXEL:
								if (reply==1)
								{
									highestPval = Pval; //global copy of max
									//apply changes
									if (chlistlen>0)
									{
										for (int q = 0; q < chlistlen; ++q)
										{
											int mm = chm[q];
											short nv = (short)newv[q];
											
											
											short[] pixels= (short[])ModelsSt.getPixels(mm+1);
											int ov = pixels[ chx[q] + modelsize*chy[q] ];
											pixels[ chx[q] + modelsize*chy[q] ] = (short)newv[q];
											
											short[] diffpx = (short[])DiffSt.getPixels(mm+1);
											int olddiffval = diffpx[chx[q] + modelsize*chy[q]];
											int newdiffval = olddiffval + nv - ov;
											diffpx[ chx[q] + modelsize*chy[q] ]=(short)(newdiffval);
											
											unsc[chm[q]] = 0;
										}
										unscwr = 0;
									}
									oldsum = newsum;
									for(int i = 0; i < use_threads; ++i)
									{
										oldsums[i] = newsums[i];
									}
									if(debug_level > 0)
									{
										Pvalhistory.print("" + looptotal);
										if(debug_level > 2)
										{
											for(int i = 0; i < use_threads; ++i)
											{
												Pvalhistory.print("\t" + newsums[i]);
											}
										}
										Pvalhistory.println("\t" + newsum);
										Pvalhistory.flush();
									}
								}
								else if( (strategy == Strategy.PIXELADJ) || (strategy == Strategy.PIX_MOVES) )  //case reply<=0
								{
									sym_factors[chm[0]] = last_factor; //restore changed sym_factors;
									
									if ( (chlistlen > 0) && (new_wraps > 0)) //dont die during warm up wrap
									{
										++unsc[chm[0]]; //must all be the same model anyways
										//for (int q=0; q<chlistlen; ++q)
										//{}
										if (unsc[chm[0]] > minfornextm &&
											( dead[chm[0]] == ModelState.CONTINUE ||
											  dead[chm[0]] == ModelState.FIXED ) )
										{
											++deadModelCount;
											System.out.println("Model: "+chm[0]+" stable for "+unsc[ chm[0] ]+" trials");
											IJ.log("Finished model " + (chm[0] + 1) ); //Frames are counted from 1
											dead[chm[0]] = ModelState.FINISHED;
											renew = 1; //request a new strategy

										}
										
										++unscwr;
										boolean byebye = false;

										if( ( target_wraps > -1) && ( wrapnum >= target_wraps) )
										{ IJ.log("ERROR: reached target for number of wraps again"); byebye = true;}
										if(deadModelCount >= target_models) //weight polishing should not be necessary
										{ IJ.log("reached target for finished models"); byebye = true;}
										if(deadModelCount == modelnum)
										{ IJ.log("all models finished"); byebye = true;}
										if(!byebye && ( unscwr > ( 3*(modelnum-deadModelCount) * minfornextm) ) )
										{
											IJ.log("no improvement for too many trials: " + unscwr );
											byebye = true;
										}

										if(byebye)
										{
											date.setTime(System.currentTimeMillis());
											IJ.log( dateFormat.format(date));
											domainloop = false;
											//break;
										}
										
									}
								}


								if (strategy == Strategy.WEIGHTADJ && !force_update)//forced weights already sit in modelweight 
								{
									if (reply>0)
									{
										for (int v=0; v<modelnum; ++v)
										{
											modelweight[v] = modelweightupd[v];
										}
									}
									else
									{
										for (int v=0; v<modelnum; ++v)
										{
											modelweightupd[v] = modelweight[v];
										}
									}
								}
								//now update the search points adjustments

								switch(strategy)
								{
									case WEIGHTADJ:
										updateWeightAdjustment();
									break;
									case PIXELADJ:
									case CONFIRMED_PIXEL:
										updatetestajustment(Models);
									break;
									default:
									break;
									//no updates required for PIXEL_EM, MODEL_EM, TOTAL_EM or PIX_MOVES
								}


								pixel_tests += chlistlen;
								chlistlen = 0;//dont reprocess in case of a forced update
								force_update = false;
								recent_save = false;
								recent_img_save = false;
							break;
							case PIXEL_EM:
								if(chlistlen!=1)
								{
									date.setTime(System.currentTimeMillis());
									System.out.println( dateFormat.format(date) );
									System.out.println("Error: chlistlen is supposed to be 1 not " + chlistlen);
									throw new IllegalStateException("Wrong chlistlen " + chlistlen);
								}
								//System.out.println("awaiting responses to Pixel_EM request");
								{
									boolean all_pixels_confirmed = true;
									for(int q = 0; q < chlistlen; ++q)
									{
										double total_mean = 0.0;
										double total_wght = 0.0;
										
										int bx = (chx[q] - homex + blocksize*(1+modelsize/blocksize) ) % blocksize;
										int by = (chy[q] - homey + blocksize*(1+modelsize/blocksize) ) % blocksize;
										int bm = (chm[q] - homem + blockdepth*(1+modelnum/blockdepth)) % blockdepth;
										int bind = bx + by * blocksize + bm * blockarea;
										double[] meanB = new double[use_threads];
										double[] wghtB = new double[use_threads];
										int[] meanS = new int[use_threads];
										int[] wghtS = new int[use_threads];
										//read the values for that expectation step
										for(int i = 0; i < use_threads; ++i)
										{
											
											int worker = threadVec.get(i).worker;
											s = reader[i].readLine();
											//System.out.println(threadVec.get(i).command);
											//System.out.println(s);
											if (s == null) // report error and return
											{
												date.setTime(System.currentTimeMillis());
												System.out.println( dateFormat.format(date) );
												System.out.println("NullString worker: " + worker + getStrategyTag() + (incrs[ind3d]) + 
												" aka: " + threadVec.get(i).command);
												IJ.log("NullString worker: " + worker + getStrategyTag() + (incrs[ind3d]) + 
												" aka: " + threadVec.get(i).command);
												throw new IllegalStateException("Unexpected reply");
											}
											
											
											Scanner sc = new Scanner(s);
											boolean scan_error = !sc.hasNext();
											boolean mean_ok = false;
											boolean wght_ok = false;
											//double meanB = 0.0;
											//int meanS = 0;
											//double wghtB = 0.0;
											//int wghtS = 0; 
											try
											{
												while(sc.hasNext() && !scan_error)
												{
													String key = sc.next();
													
													if(key.equals("MEANB"))
													{
														meanB[i] = sc.nextDouble();
														mean_ok = (meanB[i] >= 0.0);
													}
													else if(key.equals("MEANS"))
													{
														meanS[i] = sc.nextInt();	
													}
													else if(key.equals("WGHTB"))
													{
														wghtB[i] = sc.nextDouble();
														wght_ok = (wghtB[i] >= 0.0);
													}
													else if(key.equals("WGHTS"))
													{
														wghtS[i] = sc.nextInt();
													}
													else
													{
														scan_error = true;
													}
												}
											}
											catch( Exception e)
											{
												e.printStackTrace();
												scan_error = true;
											}
											
											if( scan_error || !(mean_ok && wght_ok) ) //report error and return
											{
												date.setTime(System.currentTimeMillis());
												System.out.println( dateFormat.format(date) );
												IJ.log("Unexpected reply (PIXEL_EM) from Worker:" + worker);
												IJ.log(s);
												System.out.println("received: " + s);
												System.out.println("from: " + threadVec.get(i).command);
												throw new IllegalStateException("Unexpected reply");
											}
											
											
										
										}
										//rescale shifts around 0
										int avgS = 0;
										
										for(int i = 0; i < use_threads; ++i)
										{
											avgS += meanS[i] + wghtS[i];
										}
										avgS = (int)Math.round((double)avgS/(2.0 * use_threads));
										
										for(int i = 0; i < use_threads; ++i)
										{
											meanS[i] -= avgS;
											wghtS[i] -= avgS;
										}
										//sum up the replies from the workers
										for(int i = 0; i < use_threads; ++i)
										{
											double d_mean = meanB[i]*Math.pow(2,meanS[i]);
											double d_wght = wghtB[i]*Math.pow(2,wghtS[i]);
											
											if(d_mean < 0.0)
											{
												System.out.println("Error d_mean["+i+"] of EM is negative " + d_mean);
												System.out.println(s);
												throw new IllegalStateException("raw EM expectation value was negative.");
											}
											
											if(d_wght < 0.0)
											{
												System.out.println("Error d_wght["+i+"] is negative " + d_wght);
												System.out.println(s);
												throw new IllegalStateException("weight of EM was negative.");
											}
											
											
											if( !(	Double.isInfinite(d_mean) || Double.isInfinite(d_wght) ||
													Double.isNaN(d_mean) || Double.isNaN(d_wght) 	) 
											  )
											{
												int job = threadVec.get(i).job;
												double modmin = jobVec.get(job).modelmin;
												double modmax = jobVec.get(job).modelmax;
												double dmin = jobVec.get(job).datamin;
												double dj_mean = (d_mean/d_wght + dmin - modmin)/(modmax-modmin) * modelGraylevels;
												if(dj_mean < 0.0)
												{
													System.out.println("Warning dj_mean["+i+"] of EM is negative: " + dj_mean);
													System.out.println(s + " avgS: " + avgS);
													//throw new IllegalStateException("scaled EM expectation value was negative.");
												}
												
												
												total_mean += d_wght * dj_mean;
												total_wght += d_wght;
											}
											else
											{
												System.out.println("Double overflow when processing reply from " + threadVec.get(i).command);
												System.out.println(s);
											}
										}
										//calculate to new estimate
										double nv = total_mean / total_wght;
										if(Double.isInfinite(nv) || Double.isNaN(nv))
										{
											System.out.println("Double overflow when nv = " +total_mean + "/" + total_wght);
										}
										
										int new_val = (int) ( nv + 0.5 );
										if(new_val == modelGraylevels) {--new_val;} //extremely unlikely
										if(new_val >= modelGraylevels)
										{
											System.out.println("Change: " + oldval[bind] + " -> " + new_val + " to intensity: " + nv + " is out of range");
											if(new_val > 2*modelGraylevels )
											{	throw new IllegalStateException("Averaged Pixel expectation values exceeds doubled range");}
											new_val = modelGraylevels-1;
										}
										if(new_val < 0) //no serious issue
										{
											System.out.println("Change: " + oldval[bind] + " -> " + new_val + " to intensity: " + nv + " is out of range");
											if(new_val < -modelGraylevels )
											{	throw new IllegalStateException("Averaged Pixel expectation values is too negative");}
											new_val = 0;
										}
										//System.out.println("EM_step: " + oldval[bind] + " -> " + new_val + " to intensity: " + nv);
										incrs[bind] = -100; // this value was expectation maximized
										testval[bind] = new_val;
										//System.out.println("PIXEL_EM m: " + chm[q] + " x: " + chx[q] + " y: " + chy[q] + " change: " + oldval[ind3d] + " -> " + new_val); 	
										all_pixels_confirmed &= (oldval[ind3d] == new_val);
									}
									if(all_pixels_confirmed)
									{
										strategy = Strategy.CONFIRMED_PIXEL;
										//System.out.println("switching to CONFIMED_PIXEL");
										reply = 0;
										continue;
									}
									else
									{
										strategy = Strategy.PIXELADJ;
										//System.out.println("switching to PIXELADJ");
									}
								}
								
								break;
							case MODEL_EM:
								fm = chm[0];
								lm = chm[0];
							case TOTAL_EM:
								{
									///DEBUG create backup/debug copy of current models here
									/*ImagePlus last_Models = Models.duplicate();
									last_Models.setTitle("last_Models");
									last_Models.show();
									*/
									
									int total_changes = 0;
									for(int m = fm; m <= lm; ++m)
									{
										int modelarea = modelsize * modelsize;
										double total_wght = 0.0;
										double[] wghtB = new double[use_threads];
										int[] wghtS = new int[use_threads];
								
										//read the headers for the emodels something like
										//printf("EMODEl\t%d\tPIXELS\t%d\tWGHTB\t%lf", m, modelarea, wghtB[m]);
										//printf("\tWGHTS\t%d\n", wghtS[m]); //MPFR optional
										for(int i = 0; i < use_threads; ++i)
										{
											
											int worker = threadVec.get(i).worker;
											s = reader[i].readLine();
											//System.out.print(threadVec.get(i).command+": ");
											//System.out.println(s);
											if (s == null) // report error and return
											{
												date.setTime(System.currentTimeMillis());
												System.out.println( dateFormat.format(date) );
												System.out.println("NullString worker: " + worker + getStrategyTag() + (incrs[ind3d]) +
												 " aka: " +threadVec.get(i).command);
												IJ.log("NullString worker: " + worker + getStrategyTag() + (incrs[ind3d]) +
												 " aka: " + threadVec.get(i).command);
												throw new IllegalStateException("Unexpected reply");
											}
											
											
											Scanner sc = new Scanner(s);
											boolean scan_error = !sc.hasNext();
											boolean wght_ok = false;
											int m_chk = m;
											int pix_chk = modelarea;
											try
											{
												while(sc.hasNext() && !scan_error)
												{
													String key = sc.next();
													
													if(key.equals("EMODEL"))
													{
														m_chk = sc.nextInt();
														scan_error = (m_chk != m);
														if(scan_error)
														{	System.out.println("bad EMODEL: " + m_chk + " != " + m);}
													}
													else if(key.equals("PIXELS"))
													{
														pix_chk = sc.nextInt();
														scan_error = (pix_chk != modelarea);
														if(scan_error)
														{	System.out.println("bad PIXELS: " + pix_chk + " != " + modelarea);}
													}
													else if(key.equals("WGHTB"))
													{
														wghtB[i] = sc.nextDouble();
														wght_ok = true;
													}
													else if(key.equals("WGHTS"))
													{
														wghtS[i] = sc.nextInt();
													}
													else
													{
														System.out.println("bad key: >" + key +"<");
														scan_error = true;
													}
												}
											}
											catch( Exception e)
											{
												e.printStackTrace();
												scan_error = true;
											}
											
											
											
											if( scan_error || !wght_ok ) //report error and throw
											{
												date.setTime(System.currentTimeMillis());
												System.out.println( dateFormat.format(date) );
												if(!wght_ok && !scan_error)
												{	System.out.println("bad WGHTB");}
												IJ.log("Unexpected reply (PIXEL_EM) from Worker:" + worker + " aka: " + threadVec.get(i).command);
												IJ.log(s);
												System.out.println("received: " + s);
												System.out.println("from: " + threadVec.get(i).command);
												throw new IllegalStateException("Unexpected reply");
											}
										} //end reading header lines for a model
										
										int avgS = 0;
										for(int i = 0; i < use_threads; ++i)
										{
											avgS += wghtS[i];
										}
										avgS = (int)Math.round((double)avgS/(use_threads));
										
										for(int i = 0; i < use_threads; ++i)
										{
											wghtS[i] -= avgS;
										}
										
										for(int i = 0; i < use_threads; ++i)
										{
											double d_wght = wghtB[i]*Math.pow(2,wghtS[i]);
											if( !(	Double.isInfinite(d_wght) || Double.isNaN(d_wght) ) )
											{
												total_wght += d_wght;
											}
											else
											{
												System.out.println("Double overflow when processing emodel wght from " + threadVec.get(i).command);
												System.out.println(s);
											}
										}
										System.out.println("model: " + m + " avgS: " + avgS + " weight: " + (total_wght/(use_frames)) );
										//modelweight[m] = total_wght/use_threads; //we normalize anyways but this should already be between 0 and 1
										short[] diffpx = (short[])DiffSt.getPixels(m+1);
										short[] pixels = (short[])Models.getStack().getPixels(m+1);
										for(int pix = 0; pix < modelarea; ++pix)
										{
											double total_mean = 0.0;
											double[] meanB = new double[use_threads];
											int[] meanS = new int[use_threads];
											//read the same next pixel from all different workers
											for(int i = 0; i < use_threads; ++i)
											{
												
												int worker = threadVec.get(i).worker;
												s = reader[i].readLine();
												//System.out.print(threadVec.get(i).command+": ");
												//System.out.println(s);
												if (s == null) // report error and return
												{
													date.setTime(System.currentTimeMillis());
													System.out.println( dateFormat.format(date) );
													System.out.println("NullString worker: " + worker + getStrategyTag() + (incrs[ind3d]) +
													 " aka: " + threadVec.get(i).command);
													IJ.log("NullString worker: " + worker + getStrategyTag() + (incrs[ind3d]) +
													 " aka: " + threadVec.get(i).command);
													throw new IllegalStateException("Unexpected reply");
												}
												
												
												Scanner sc = new Scanner(s);
												boolean scan_error = !sc.hasNext();
												try
												{
													meanB[i] = sc.nextDouble();
													if(sc.hasNext())
													{	meanS[i] = sc.nextInt();}
													meanS[i] -= avgS;  
												}
												catch( Exception e)
												{
													e.printStackTrace();
													scan_error = true;
												}
												
												if( scan_error  ) //report error and return
												{
													date.setTime(System.currentTimeMillis());
													System.out.println( dateFormat.format(date) );
													IJ.log("Failed to read Model: " + m + " Pixel: "  + pix + " from: " + threadVec.get(i).command);
													IJ.log(s);
													System.out.println("received: " + s);
													System.out.println("from: " + threadVec.get(i).command);
													throw new IllegalStateException("Unexpected reply" + s + " aka " + threadVec.get(i).command);
												}
											}// end reading the same pixel from all workers
											//only further process those inside the hexagon
											//int q = pix % modelsize;
											//int r = pix / modelsize;
											
											if(hex_mask[pix])//(iswithinHex(q,r,modelsize/2))
											{
												//sum up the pixel_values from the workers
												for(int i = 0; i < use_threads; ++i)
												{
													double d_mean = meanB[i]*Math.pow(2,meanS[i]);
													double d_wght = wghtB[i]*Math.pow(2,wghtS[i]);
													if( !(	Double.isInfinite(d_mean) || Double.isNaN(d_mean) ) )
													{
														int job = threadVec.get(i).job;
														double modmin = jobVec.get(job).modelmin;
														double modmax = jobVec.get(job).modelmax;
														double dmin = jobVec.get(job).datamin;
														double dj_mean = (d_mean/d_wght + dmin - modmin)/(modmax-modmin) * modelGraylevels;
														total_mean += d_wght * dj_mean;
													}
													else
													{
														date.setTime(System.currentTimeMillis());
														System.out.println( dateFormat.format(date) );
														System.out.println("Double overflow when processing reply from " + threadVec.get(i).command);
														System.out.println(s);
														throw new IllegalStateException("Math Error");
													}
												}
												
												double nv = total_mean / total_wght;
												if(Double.isInfinite(nv) || Double.isNaN(nv))
												{
													System.out.println("Double overflow when nv = " + total_mean + "/" + total_wght);
												}
												int new_val = (int) ( nv + 0.5 );
												if(new_val == modelGraylevels) {--new_val;} //extremely unlikely
												if(new_val >= modelGraylevels)
												{
													System.out.println("Change: " + pixels[pix] + " -> " + new_val + " to intensity: " + nv + " is out of range");
													if(new_val > 2 * modelGraylevels)
													{	throw new IllegalStateException("Averaged Pixel expectation value exceeds doubled range");}
													new_val = modelGraylevels-1;	
												}
												if(new_val < 0)
												{	new_val = 0;}
												total_changes += (int)Math.abs((short)(new_val - pixels[pix]));
												diffpx[pix] += (short)(new_val - pixels[pix]);
												pixels[pix] = (short)new_val;
											}
										
										}// end for pix
									}//end for m
									//normalize(modelweight);
									if(total_changes > 0)
									{
										for(int i = 0; i < use_threads; ++i)
										{
											//pwos[i].println("updateWeights()");
											//pwos[i].flush();
											//sendWeights(pwos[i], os[i], modelweight);
											pwos[i].println("updateModels("+fm+","+lm+")");
											pwos[i].flush();
											sendModels(pwos[i], os[i], Models);
											pwos[i].println("Run(1)");
											pwos[i].flush();
										}
										strategy = Strategy.PIXELADJ;
										force_update = true;
										ongoing_Adj = -1;
										//force_wrap = true;
										continue; //with force_update loop
									}
									else
									{
										strategy = Strategy.PIXELADJ;
										ongoing_Adj = -1;
										do_full_EM = false;
										IJ.log("No further changes suggested by EM alogrithm, cancelling further wraps");
										wraps_completed = true;
									}
									break;
								}	
							case POLL_HISTO:
							{
								++histonum;
								IJ.showStatus("Polling Histo: " + histonum);
								
								double[][] pixels = new double[ jobVec.size() ][];
								 
								
								for(int i = 0; i < use_threads; ++i)
								{
									final int job = threadVec.get(i).job;
									final int ptwidth = 1+jobVec.get(job).datamax-jobVec.get(job).datamin;
									final int ptableSize = ptwidth * modelGraylevels;
									
									if( pixels[job] == null)
									{
										//IJ.log("allocated pixels for job: " + job);
										pixels[job] = new double[ptableSize];
									} 
									
									
									int worker = threadVec.get(i).worker;
									s = reader[i].readLine();
									//System.out.print(threadVec.get(i).command+": ");
									//System.out.println(s);
									if (s == null) // report error and return
									{
										date.setTime(System.currentTimeMillis());
										System.out.println( dateFormat.format(date) );
										System.out.println("NullString worker: " + worker + getStrategyTag() + (incrs[ind3d]) +
										 " aka: " +threadVec.get(i).command);
										IJ.log("NullString worker: " + worker + getStrategyTag() + (incrs[ind3d]) +
										 " aka: " + threadVec.get(i).command);
										throw new IllegalStateException("Unexpected reply");
									}
									
									
									Scanner sc = new Scanner(s);
									boolean scan_error = !sc.hasNext();
									boolean histo_ok = false;
									int num_points = -1;
									boolean base_ok = false;
									boolean shift_ok = false;
									
									try
									{
										while(sc.hasNext() && !scan_error)
										{
											String key = sc.next();
											
											if(key.equals("HISTO"))
											{
												histo_ok = true;
											}
											else if(key.equals("POINTS"))
											{
												num_points = sc.nextInt();
												scan_error = (num_points != ptableSize);
												if(scan_error)
												{	System.out.println("wrong number of POINTS: " + num_points + " != " +  ptableSize);}
											}
											else if(key.equals("BASE"))
											{
												base_ok = true;
											}
											else if(key.equals("SHIFT"))
											{
												shift_ok = true;
											}
											else
											{
												System.out.println("bad key: >" + key +"<");
												scan_error = true;
											}	
										}
										if(!(histo_ok && base_ok))
										{
											scan_error = true;	
										}
									}
									catch( Exception e)
									{
										e.printStackTrace();
										scan_error = true;
									}
									
									if( scan_error ) //report error and throw
									{
										date.setTime(System.currentTimeMillis());
										System.out.println( dateFormat.format(date) );
										IJ.log("Unexpected reply (updateHisto) from Worker:" + worker + " aka: " + threadVec.get(i).command);
										IJ.log(s);
										System.out.println("received: " + s);
										System.out.println("from: " + threadVec.get(i).command);
										throw new IllegalStateException("Unexpected reply" + s + " aka: " + threadVec.get(i).command);
									}
									try
									{
										for(int ind = 0; ind < ptableSize; ++ind)
										{
											double base = 0;
											int shift = 0;
											s = reader[i].readLine();
											Scanner ls = new Scanner(s);
											base = ls.nextDouble();
											if(Double.isNaN(base) || Double.isInfinite(base))
											{
												throw new NumberFormatException(s);
											}
											if(shift_ok)
											{
												shift = ls.nextInt();
												base *= Math.pow(2,shift);
												if(Double.isNaN(base) || Double.isInfinite(base))
												{
													throw new NumberFormatException(s);
												} 
											}
											
											pixels[job][ind] += base; // //random.nextDouble() does not help against all zeros
										}
									}
									catch( Exception e)
									{
										System.out.println(threadVec.get(i).command);
										System.out.println("failed to parse Values from: " + s );
										//make sure the entire histogram is purged
										while(reader[i].ready())
										{
											while(reader[i].ready())
											{
												reader[i].readLine();
											}
											Thread.sleep(25); 
										}		
									}
								}
								
								for( int job = 0; job < jobVec.size(); ++job)
								{
									
									double pvalfloor = - 3.0;
									final int ptwidth = 1+jobVec.get(job).datamax-jobVec.get(job).datamin;
									final int ptableSize = ptwidth * modelGraylevels;
									String nptitle = "job" + job + "_Histo.tif";
									ImagePlus NewHisto = WindowManager.getImage(nptitle);
									if(NewHisto == null)
									{
										//IJ.log("created " + nptitle);
										NewHisto = NewImage.createFloatImage( nptitle , ptwidth , modelGraylevels , 1, NewImage.FILL_BLACK);
										ImageStack hst = NewHisto.getStack();
										hst.setSliceLabel("wrap_"+wrapnum,1);
										NewHisto.setStack(hst);
										NewHisto.show();
									}
									ImageStack NewHistoStack = NewHisto.getStack();
									while(NewHistoStack.getSize() < histonum)
									{
										float[] tmp_pixels = new float[ptableSize];
										//IJ.log("adding slice to " + nptitle);
										NewHistoStack.addSlice("wrap_" + wrapnum, tmp_pixels);
									}
									NewHisto.setStack(NewHistoStack);
									float[] fpixels = (float[])NewHistoStack.getPixels(histonum);
									
									double[] area = new double[modelGraylevels];
									for(int ind=0; ind < ptableSize; ++ind)
									{
										int y = ind/ptwidth;
										area[y] += (pixels[job][ind]);
									}
									for(int y = 0; y < modelGraylevels; ++y)
									{
										final int i1 = y * ptwidth;
										final int i2 = i1 + ptwidth; 
										if(area[y] > 0)
										{
											for(int ind = i1; ind < i2; ++ind)
											{
												fpixels[ind] = (float)Math.log(pixels[job][ind]/area[y]);	
											}
										}
										else
										{
											for(int ind = i1; ind < i2; ++ind)
											{
												fpixels[ind] = (float)(2*pvalfloor);
											}
										}	
									}
									for(int ind = 0; ind < ptableSize; ++ind)
									{
										if(fpixels[ind] < pvalfloor)
										{
											fpixels[ind] = (float) (pvalfloor*(1.0d+Math.tanh((fpixels[ind]-pvalfloor)/pvalfloor)));
										}
									}
									NewHisto.updateAndRepaintWindow();
								}
								//TODO normalize fpixels to same avg as Ptable and decide wheter to employ it or not
								
							}
							break;
							default:
								date.setTime(System.currentTimeMillis());
								System.out.println( dateFormat.format(date) );
								throw new IllegalStateException("Undefined STRATEGY when processing feedback from workers");
									
						}
						//delay esc action until there is no other output printed
						boolean loopclear = ( loopcounter <= loopmax - report_num + 1);
						loopclear |= (loopmax <= report_num); //make sure esc is accepted in pathological cases
						final boolean pause_now = ( pause_at_next_wrap && (wrapnum == pause_at_wrap) );
						if( pause_now || (IJ.escapePressed() && loopclear) )
						{	
							//TODO turn this duplicate code a function
							if( !( (debug_level > 1) && exclusive_autosave ) )
							{
								System.out.println(master_path);
								saveOptions(options, false);
								System.out.println("saved " + options.getName());
								saveWeights(master_path);
								if( can_do_stats )
								{	
									writeSignificance(master_path, Models);
									System.out.println("saved " + "Significance.txt");
								}
								recent_save = true;
								saveAllImgs( master_path );
								Process cpPval = rt.exec("cp -f " + master_path + "Pval.txt " + master_path + "Pval.bak"); 
								cpPval.waitFor();
								cpPval.destroy();
								System.out.println("Pval.txt -> Pval.bak");
								//Modeltitle = Models.getTitle();
								recent_img_save = true;
							}
							
							IJ.resetEscape();
							System.out.println( getClass().getSimpleName() + " halted");
							IJ.open(master_path + "Weights.txt");
							IJ.selectWindow("Weights.txt");
							NonBlockingGenericDialog nbgd = new NonBlockingGenericDialog( getClass().getSimpleName() + " halted");
							nbgd.addMessage("Optimzation is paused loopcount = " + looptotal);
							nbgd.addMessage("You may now edit" + master_path + "Weights.txt");
							nbgd.addCheckbox("reload weights", false);
							nbgd.addCheckbox("reinitialize models (BEWARE!)",false);
							nbgd.addCheckbox("resend models", false);
							nbgd.addCheckbox("resend weights", false);
							nbgd.addMessage("Changes in the following settings may require resending and/or reloading");
							nbgd.addCheckbox("(re-)normalize model intensities", false);
							nbgd.addCheckbox("noisify models", noisify_at_restart);
							nbgd.addNumericField("noise area diameter:", modelactiveinit, 0);
							nbgd.addMessage("The settings below only affect the master");
							nbgd.addMessage("The new pattern blocksize and blockdepth will only be applied when the nest wrap starts");
							
							nbgd.addChoice("scanning pattern (next wrap)", pattern, scan_pattern);
							nbgd.addNumericField("blocksize:", new_blocksize, 0);
							nbgd.addNumericField("blockdepth:", new_blockdepth, 0);
							
							nbgd.addNumericField("target models:", target_models, 0);
							nbgd.addNumericField("target wraps:", target_wraps, 0);
							nbgd.addCheckbox("fixed_lattice_weight", fixed_lattice_weight);
							nbgd.addCheckbox("frozen_lattice", frozen_lattice);
							nbgd.addCheckbox("use weights (next wrap)", next_use_weights);
							nbgd.addCheckbox("const weights (next wrap)", const_weights); //anyway only effective after a wrap
							nbgd.addCheckbox("use max (next wrap)", next_use_max);
							
							nbgd.addCheckbox("alternating scan orientation", alternating_scan_direction);
							nbgd.addCheckbox("flip block", flip_block);
							nbgd.addCheckbox("blockdepth priority", blockdepth_priority);
							nbgd.addCheckbox("reverse block scan", reverse_block_depth);
							
							nbgd.addCheckbox("correlated optimization", correlated_optimization);
							nbgd.addCheckbox("quick optimization", quick_optimization);
							nbgd.addCheckbox("rough optimization", rough_optimization);
							nbgd.addCheckbox("enable quad fits", enable_quad_fits);
							nbgd.addCheckbox("black initalization", black_init);
							nbgd.addCheckbox("exclusive autosave", exclusive_autosave);
							nbgd.addCheckbox("blobs (very expensive!)", blobs);
							//This will only work in normal reconstructions with PIXEL updates and actual weights
							nbgd.addCheckbox("pause again right before next wrap",  pause_at_next_wrap );
							nbgd.addCheckbox("clone models (after the run)", clone_models);
							nbgd.addCheckbox("but dont clone first model", skip_first_clone);
							nbgd.addCheckbox("interactive cloning session",interactive_cloning_session);
							nbgd.addNumericField("rlative weight for clones", clone_weight,3);
							
							nbgd.addMessage("debug level " + debug_level);
							nbgd.addMessage("Click Ok to resume");

							nbgd.showDialog();
							if(WindowManager.getWindow("Weights.txt") != null)
							{
								IJ.selectWindow("Weights.txt");
								IJ.run("Close");
							}
							if(nbgd.wasCanceled())
							{
								domainloop = false;
								clone_models = false; //only do that at a clean finish
								interactive_cloning_session = false;
							}
							else
							{
								boolean reload_weights = nbgd.getNextBoolean();
								boolean reinit_models = nbgd.getNextBoolean();
								boolean resend_models = nbgd.getNextBoolean();
								boolean resend_weights = nbgd.getNextBoolean();
								
								normalize_model_intensities = nbgd.getNextBoolean();
								noisify_at_restart = nbgd.getNextBoolean();
								
								modelactiveinit = (int)nbgd.getNextNumber();
								modelactive = modelactiveinit; //DUPLICATE CODE HERE
								if(modelactive <= -60)
								{	modelactive = - (60 * modelsize) / modelactive;}
								else if(modelactive < 0)
								{	modelactive = -bondlength * modelactive;}
								if(modelactive > modelsize)
								{	modelactive = modelsize;}
								
								scan_pattern = nbgd.getNextChoice();
								new_blocksize = (int)nbgd.getNextNumber();
								new_blockdepth = (int)nbgd.getNextNumber();
								
								target_models = (int)nbgd.getNextNumber();
								target_wraps = (int)nbgd.getNextNumber();
								fixed_lattice_weight = nbgd.getNextBoolean();
								frozen_lattice = nbgd.getNextBoolean();
								next_use_weights = nbgd.getNextBoolean();
								const_weights = nbgd.getNextBoolean();
								next_use_max = nbgd.getNextBoolean();
								
								alternating_scan_direction = nbgd.getNextBoolean();
								flip_block = nbgd.getNextBoolean();
								blockdepth_priority = nbgd.getNextBoolean();
								reverse_block_depth = nbgd.getNextBoolean();
								
								correlated_optimization = nbgd.getNextBoolean();
								quick_optimization = nbgd.getNextBoolean();
								rough_optimization = nbgd.getNextBoolean();
								enable_quad_fits = nbgd.getNextBoolean();
								boolean next_black_init = nbgd.getNextBoolean();
								exclusive_autosave = nbgd.getNextBoolean();
								blobs = nbgd.getNextBoolean();
								pause_at_next_wrap = nbgd.getNextBoolean();
								clone_models = nbgd.getNextBoolean();
								skip_first_clone = nbgd.getNextBoolean();
								interactive_cloning_session = nbgd.getNextBoolean();
								clone_weight = nbgd.getNextNumber();
								if(pause_at_next_wrap)
								{	pause_at_wrap = wrapnum + 1; }
								if(new_blocksize < 0)
								{	new_blocksize = -new_blocksize;}
								if(new_blockdepth < 0)
								{	new_blockdepth = modelnum;}
								if(new_blocksize > modelsize)
								{	new_blocksize = modelsize;}
								if( new_blockdepth > modelnum)
								{	new_blockdepth = modelnum;}
								
								
								if(reload_weights)
								{
									while ( !initialize_weights(master_path, Weighttitle))
									{
										IJ.error("Please fix " + master_path + Weighttitle);
									}
									if(!resend_weights)
									{
										IJ.showMessage("New weights will be broadcasted to workers now!");
									}
									resend_weights = true;
								}
								if(reinit_models)
								{
									initialize_models( Models, Diffimg, modelGraylevels, modelmean0);
									ModelsSt = Models.getStack();
									Models.updateAndDraw();
									Diffimg.updateAndDraw();
									if(!resend_models)
									{
										IJ.showMessage("New models will be broadcasted to workers now!");
									}
									resend_models = true;
								}
								else if(normalize_model_intensities)
								{
									IJ.showMessage("re-normalization Models would require re-initalization!");
									normalize_model_intensities = false;
								}
								if(resend_models)
								{	System.out.println("resending models"); }
								for(int i = 0; i < streamnum; ++i)
								{
									boolean finalize = false;
									if(next_black_init != black_init)
									{
										pwos[i].println("black_init(" + (black_init?1:0)  + ")");
										pwos[i].flush();
									}
									
									if(resend_weights)
									{
										pwos[i].println("updateWeights()");
										pwos[i].flush();
										sendWeights( pwos[i], os[i], modelweight);
										finalize = true;
										//force_update = true; //they require an acknowledgement as they are like an regular update
									}
									if(resend_models)
									{
										pwos[i].println("updateModels()");
										pwos[i].flush();
										sendModels(pwos[i], os[i], Models);
										finalize = true;
										//force_update = false; //a forced reinitialization never requires acknowledgements
									}
									if(finalize) //some data was actually updated
									{
										pwos[i].println("Run(1)");
										pwos[i].flush();
										force_update = true;
									}
								}
								black_init = next_black_init;
								System.out.println("resuming " + getClass().getSimpleName() );
							}	
						}
					} while (force_update && domainloop);
					loopcounter++;
					looptotal++;
					sync_images( new ImagePlus[]{Models, Diffimg, Updates, Activity} );
				} while (domainloop);

				for(int i=0; i < streamnum; ++i)
				{
					if(pending_keep)
					{	pwos[i].println("Keep(-1)");}
					pwos[i].println("Run(0)");
					pwos[i].flush();
				}
				pending_keep = false;
				for (int i=0; i<use_threads; ++i)
				{
					//System.out.println("Shutting down worker " + i);
					p[i].waitFor();
					p[i].destroy();
				}
				{
					Process ps4 = rt.exec( "rm " + nuke.getPath() );
					ps4.waitFor();
					ps4.destroy();
				}
			} 
			catch(Exception e)
			{
				e.printStackTrace();
				for (int i=0; i<streamnum; ++i)
				{
					if(pending_keep)
					{
						pwos[i].println("Keep(1)");
					}
					if(can_write[i])
					{	
						int worker = threadVec.get(i).worker;
						System.out.println("Requesting writeout for cachereport" + worker + ".txt");
						pwos[i].println("CacheReport(-1)");
					}
					
					pwos[i].println("Run(0)");
					pwos[i].flush();
				}
				boolean clean_shut_down = true;
				for (int i=0; i<use_threads; ++i)
				{
					
					try
					{	p[i].waitFor();}
					catch (Exception e2)
					{
						int worker = threadVec.get(i).worker;
						clean_shut_down = false;
						e2.printStackTrace();
						//actually happens if an exception occurs while writing job files
						System.out.println("Error while shutting down worker" + worker);
					}
					p[i].destroy();
				}
				
				if(clean_shut_down)
				{
					try
					{
						if( nuke.isFile())
						{
							Process ps5 = rt.exec( "rm " + nuke.getPath() );
							ps5.waitFor();
							ps5.destroy();
						}
					}
					catch (Exception e2nd)
					{
						System.out.println("WARNING: failed to remove temporary nuke script");
						e2nd.printStackTrace();
					}
				}
				
				IJ.showMessage("An error has occured!");
			}
			if(clone_models)
			{
				modelManipulation(false);
			}
			if(!batch_run)
			{
				date.setTime(System.currentTimeMillis());
				IJ.log( dateFormat.format(date));
			}
			if( clone_models || !( (debug_level > 1) && exclusive_autosave ) )
			{
				if(clone_models || !recent_save)
				{
					System.out.println( master_path);
					if( !initbenchmark)
					{
						saveWeights(master_path);
						saveOptions(options, false);
						System.out.println("saved " + options.getName());
					
						if(can_do_stats)
						{
							writeSignificance(master_path, Models);
							System.out.println("saved " + "Significance.txt");
						}
					}
				}
				if( (clone_models || !recent_img_save) && !initbenchmark)
				{			
					saveAllImgs( master_path );
					try
					{
						Process cpPval = rt.exec("cp -f " + master_path + "Pval.txt " + master_path + "Pval.bak"); 
						cpPval.waitFor();
						cpPval.destroy();
						System.out.println("Pval.txt -> Pval.bak");
					}
					catch(Exception e)
					{
						e.printStackTrace();
						//actually no serious problem
					}	
				}
			}
			saveLog(master_path);
			System.out.println("saved " + "Log.txt");
			if(!batch_run)
			{
				IJ.open(master_path + "Log.txt");
			}
			if(WindowManager.getWindow("Log") != null)
			{
				IJ.selectWindow("Log");
				IJ.run("Close");
			}
			
			if(batch_run) 
			{	//close all modified and hopefully saved images
				Models.close();
				Histimg.close();
				Diffimg.close();
				if(!Modus.equals(Modii[2]))
				{
					for(int i = 0; i < jobVec.size(); ++i)
					{	
						jobVec.get(i).Ptable.close();
						if(!keep_Idata)
						{	jobVec.get(i).Idata.close();}
					}
				}
				if(Activity != null) Activity.close();
				if(Updates != null) Updates.close();
				if(WindowManager.getWindow("batchrun.txt") == null)
				{
					statusdisplay.append("batchrun.txt was closed, skipping all further runs");
					break; //cancel batch mode
				}
			}
		} //endfor batch_id   
		if(batch_run)
		{
			statusdisplay.append("batchrun completed");
		}
		IJ.showStatus( getClass().getSimpleName() + " has finished");
		IJ.showProgress(1.0);
		System.out.println( getClass().getSimpleName() + " has finished");
		mutex = 0;
    }

	void saveImg(ImagePlus imgp, String path)
    {
		if( imgp!=null )
		{
			if(imgp.getStackSize() > 1)
			{	new FileSaver(imgp).saveAsTiffStack( path + imgp.getTitle() );}
			else
			{	new FileSaver(imgp).saveAsTiff( path + imgp.getTitle() );}
			System.out.println("saved " + imgp.getTitle() );
		}
	}
	
	void saveAllImgs(String path)
	{
		saveImg( Models, path);
		saveImg( Diffimg, path);
		saveImg( Activity, path);
		saveImg( Updates, path);
		saveImg( Histimg, path);
		saveImg( beamProfile, path);
	}
	
	
	int[] nextModelSequence(int first, int last, int num_clones)
	{
		
		System.out.println("model cloning(" + modelnum + "), first: " + first + " last: " + last + " num_clones: " + num_clones);
		if( (num_clones < 0) || (first > last) )
		{
			System.out.println("Error invalid arguments in nextModelSequence first: " +
				first + ", last: " + last + ", num_clones: " + num_clones);
			return null;	
		}
		int nxt_modelnum = modelnum + (last-first+1)*(num_clones-1);
		
		/*
		switch(num_clones)
		{
			case 0: //deletion
			nxt_modelnum -= (last-first+1);
			break;
			case 1:
			break;
			default:
			nxt_modelnum += (last-first+1)*(num_clones-1);
		}
		*/ 
		int[] nxt_seq = new int[nxt_modelnum];
		
		int mdl2 = 0;
		for(int mdl = 1; mdl <= modelnum;++mdl)
		{
			if(mdl < first || mdl > last) //regular copy
			{	nxt_seq[mdl2++]=mdl;}
			else 
			{
				
				for(int ncls = 0; ncls < num_clones; ++ncls) //multiple copy (including none)
				{
					nxt_seq[mdl2++]=mdl;
				}
				
			}
		}
		if(num_clones == 1) //swap 
		{
			nxt_seq[first-1] = last;
			nxt_seq[last-1] = first;
		}
		System.out.print("next model sequenxe("+nxt_modelnum+"): {");
		for( int m = 0; m < nxt_modelnum; ++m)
		{
			System.out.print(" " + nxt_seq[m]);
		}
		System.out.println(" }");
		return nxt_seq;
	}
	
	
	
	ImagePlus expandImg(ImagePlus imp, int [] nxt_seq) //one based indices
	{
		
		
		ImagePlus orig = imp; //backup copy
		ImageStack origSt = orig.getStack();
		int w = orig.getWidth();
		int h = orig.getHeight();
		int slices = orig.getStackSize();
		int newmodels = nxt_seq.length - modelnum;
		
		if(imp.getType() == ImagePlus.GRAY16) //Models, Diff, Activity or Updates
		{
			int area = w*h;
			imp = NewImage.createShortImage(orig.getTitle(), w, h ,nxt_seq.length, NewImage.FILL_BLACK);
			ImageStack impSt = imp.getStack();
			int slo = 0;
			
			for(int sl = 1; sl <= nxt_seq.length; ++sl)
			{
				
				slo = nxt_seq[sl-1];
				short[] imppix2  = (short[]) impSt.getPixels(sl);
				short[] origpix = (short[])origSt.getPixels(slo);
				for(int ind = 0; ind < area; ++ind)
				{	imppix2[ind] = origpix[ind];}	
			}
			imp.setStack(impSt);
		}
		else if (imp.getType() == ImagePlus.GRAY32) //must be History
		{
			int mdlsz = w/modelnum;
			int nw = (modelnum + newmodels)*mdlsz;
			//int left_edge = (skip_first_clone && (modelnum > 1) ) ? modelsize+1: 0;
			imp = NewImage.createFloatImage(orig.getTitle(), nw , h , slices , NewImage.FILL_BLACK);
			ImageStack impSt = imp.getStack();
			for(int sl = 1; sl <= slices; ++sl)
			{
				impSt.setSliceLabel(origSt.getSliceLabel(sl),sl);
				float[] origpix = (float[])origSt.getPixels(sl);
				float[] imppix  = (float[]) impSt.getPixels(sl);
				for(int mdl = 1; mdl <= nxt_seq.length; ++mdl)
				{
					int mdl2 = (mdl-1);
					int mdl1 = nxt_seq[mdl2]-1;
					
					int left1 = mdl1*mdlsz; 
					int left2 = mdl2*mdlsz;
					
					for(int y=0; y < h; ++y)
					{
						int pos1 = y*w+left1;
						int pos2 = y*nw+left2;
						for(int x = 0; x < mdlsz; ++x)
						{	
							imppix[pos2++] = origpix[pos1++];
						}	
					} 
				} 	
			}
			imp.setStack(impSt);	
		}
		orig.close(); //comment this to keep the previous version
		imp.show();
		return imp;
	}
    
    void expandWeights( int first, int last, int num_clones )
    {
		if (first > last)
		{	
			System.out.println("Error in expandWeights first:" + first + " > last:" + last);
			return;
		}
		if ( num_clones < 0 )
		{	
			System.out.println("Error in expandWeights num_clones: " + num_clones);
			return;
		}
		boolean cp_rm = num_clones > 0;
		int newmodels = (last-first + 1) * (num_clones-1);
		int nxt_modelnum = modelnum+newmodels;
		
		
		double[] nxt_modelweight = new double[nxt_modelnum];
		ModelState[] nxt_dead = new ModelState[nxt_modelnum];
		int mdl1 = 0;
		int mdl2 = 0;
		double cl_weight = (num_clones > 1)? clone_weight/(num_clones-1) : 0.0;
		for(int mdl=1; mdl <= modelnum; ++mdl)
		{
			if((!cp_rm) && (mdl>=first) && (mdl<=last))
			//simply skip originals if !cp_rm aka removing
			{	
				++mdl1;
				continue;
			}
			int copies = 0;
			do
			{
				if(!((mdl>=first) && (mdl<=last)))//outside selected range regular copy
				{
					nxt_modelweight[mdl2] = modelweight[mdl1];
					nxt_dead[mdl2] = dead[mdl1];
				}
				else if(cp_rm)
				{
					nxt_dead[mdl2] = ModelState.CONTINUE; //cloned models will always be active
					nxt_modelweight[mdl2] = ((copies==0) ? (1.0-cl_weight) : cl_weight)*modelweight[mdl1];
				}
				else
				{
					System.out.println("Error unexpected cp_rm: " + cp_rm + " encountered at mdl:" + mdl);
				}
				++mdl2;
			}
			while( cp_rm && (mdl>=first) && (mdl<=last) && (++copies < num_clones) ); //2 children per ancestor
			++mdl1;	
		}
		if(num_clones == 1)
		{
			nxt_modelweight[first-1] = modelweight[last-1];
			nxt_modelweight[last-1] = modelweight[first-1];
			nxt_dead[first-1] = dead[last-1]; //no special state for swapped models
			nxt_dead[last-1] = dead[first-1];
		}
		
		dead = nxt_dead;
		modelnum = nxt_modelnum;
		normalizeWeights( nxt_modelweight );
		modelweight = nxt_modelweight;
		modelweightupd = new double[modelnum];				
		unsc = new int[modelnum];
		cases = new int[modelnum];
		sym_factors = new double[modelnum];
		
		final int num_blocks = (modelsize/blocksize) * (modelsize/blocksize);
		permutation = new int[modelnum][num_blocks];
		init_perm(permutation);
		Diff_Matcher matcher = null;
		if(sym_scaling)
		{
			matcher = new Diff_Matcher(Models, 1, modelnum, bondlength, (median_filtering?modelGraylevels:0) );
		}
		if(matcher == null)
		{
			Arrays.fill(sym_factors,1.0);
		}
		else
		{	
			for(int m = 0; m < modelnum; ++m)
			{
				sym_factors[m] = median_filtering ? matcher.summed_matches[m] : matcher.inv_avg_matches[m];	
			}	
		}
	}
    
    
    
    
    void send16bitint(java.io.OutputStream ostr, short val)
    throws java.io.IOException  //low bits first.
    {
        buffer.clear();
        buffer.putShort(0, (short)val);
        ostr.write( bytes, 0, 2);
    }

    void senddouble(java.io.OutputStream ostr, double val)
    throws java.io.IOException
    {
        buffer.clear();
        buffer.putDouble(0, (double)val);
        ostr.write( bytes, 0, 8 );
    }
	//This one does not work. some issue with byte buffer
    double getdouble(java.io.InputStream istr)
    throws java.io.IOException
    {
        buffer.clear();
        istr.read( bytes, 0, 8 );
        return buffer.getDouble();
    }

	void sendWeights(java.io.PrintWriter pwos, java.io.OutputStream ostr, double[] weights)
	throws java.io.IOException
	{
		
		double weightSum = 0.0;
		pwos.println("BeginBinary(" + (8 * modelnum) + ")");
		pwos.flush();
		for (int j = 0; j < modelnum; ++j)
		{
			double scaled_weight = weights[j] * sym_factors[j];
			
			senddouble(ostr,scaled_weight);
			weightSum += scaled_weight;
		}
		ostr.flush();
		pwos.println("EndBinary(" + weightSum + ")");
		pwos.flush();
	}


	void sendModels(java.io.PrintWriter pwos, java.io.OutputStream ostr, ImagePlus Models)
	throws java.io.IOException
	{
		//send model images
		ImageStack ModelsSt = Models.getStack();
		int chsum = 0;
		int lastmodel = ModelsSt.getSize(); // 1 based
		int modelarea = Models.getWidth() * Models.getHeight();
		pwos.println("BeginBinary(" + (2 * lastmodel * modelarea) + ")");
		pwos.flush();
		for (int v=1; v<= lastmodel; ++v)
		{
			short[] pixels = (short[])ModelsSt.getPixels(v);
			for (int j=0; j<pixels.length; ++j)
			{
				send16bitint(ostr,pixels[j]);
				chsum +=(int)pixels[j] & 0xffff;
			}
			ostr.flush();
		}

		pwos.println("EndBinary(" + chsum + ")");
		pwos.flush();
	}


    double normalizeWeights(double[] vec)
	{
		double sum = 0.0;
		double lattice_weight = 0.0;
		int len = vec.length;

		for (int j = 0; j < len ; ++j)
		{
			if ( dead[j] != ModelState.LATTICE && dead[j] != ModelState.FIXED )
			{
				if(vec[j] <= 0)
				{
					vec[j] = -vec[j];
					IJ.log("WARNING, weight of model" + (j+1) + " was -"+ vec[j]);
				}
				sum += vec[j];
			}
			else
			{	lattice_weight += vec[j];}
		}

		if (sum != 0.0)
		{
			sum = (1.0-lattice_weight)/sum;
			for (int j = 0; j < len; ++j)
			{
				if (dead[j] != ModelState.LATTICE && dead[j] != ModelState.FIXED  )
				{	vec[j] *= sum;}
			}
		}
		return 1.0/sum;
	}

	double normalize(double[] vec)
	{
		double sum = 0.0;
		int len = vec.length;
		for (int j = 0; j < len ; ++j)
		{
			vec[j] = Math.abs(vec[j]);
			sum += vec[j];
		}

		if (sum != 0.0)
		{
			double invsum = 1.0/sum;
			for (int j = 0; j < len; ++j)
			{
				vec[j] *= invsum;
				
			}
		}
		return sum;
	}

	void shuffle(int[] vec)
	{
		for(int i = 1; i<vec.length; ++i)
		{
			int j = random.nextInt(i+1);
			if(i != j)
			{
				int tmp = vec[j];
				vec[j] = vec[i];
				vec[i] = tmp;
			}
		}	
	}
	
	void shuffle(int[][] vec)
	{
		for(int i = 0; i<vec.length; ++i)
		{	shuffle(vec[i]);}
	}
	
	void init_perm(int[] vec)
	{
		for(int i = 0; i<vec.length; ++i)
		{	vec[i] = i;}
		shuffle(vec);
	}
	
	void init_perm(int[][] vec)
	{
		for(int i = 0; i<vec.length; ++i)
		{	init_perm(vec[i]);}
	}
	
	//extracts an shuffled index and suffles vec
	//i should be incremented between successive calls
	int get_and_shuffle(int i,int[] vec)
	{
		i = i % vec.length;
		int j = random.nextInt(i+1);
		if(i != j)
		{
			int tmp = vec[j];
			vec[j] = vec[i];
			vec[i] = tmp;
		}
		return vec[j];
	}
	

    void initialize_models(ImagePlus Models, ImagePlus Diffimg, int modelGraylevels, double modelmean)
    {
        //IJ.log("initializing models: graylevels = " + modelGraylevels + "modelmean = " + modelmean);
        int checksum = 0;
        //The lattice has always to be the first model
        int j = 1;
        int active_radius = modelactive/2;
        String Mtitle = Models.getTitle();
        ImageStack ModelsSt = Models.getStack();
        ImageStack DiffSt = Diffimg.getStack();

        if (new_models)
        {
			short[] pixels = (short[])ModelsSt.getPixels(j);
			for(int i = 0; i < pixels.length; ++i)
			{
				//int q = i % modelsize;
				//int r = i / modelsize;
				if(hex_mask[i])// (iswithinHex(q,r,modelsize/2))
				{
					double noise = -1.0;
					while (noise < 1.0)
					{
						noise = modelmean + modelmean*area_noise_level*random.nextGaussian();
					}
					pixels[i] = (short)noise;
				}
			}
		}
		if( (dead[j-1] == ModelState.LATTICE) )
		{
			if (frozen_lattice)
			{
				if ( fixed_lattice_weight)
				{	dead[j-1] = ModelState.LATTICE;}
				else
				{	dead[j-1] = ModelState.FINISHED;}
			}
			else // !frozen_lattice
			{
				--deadModelCount;
				if ( fixed_lattice_weight)
				{	dead[j-1] = ModelState.FIXED;}
				else
				{	
					dead[j-1] = ModelState.CONTINUE;
				}
			}
			
		}
        short[] lattice = (short[])ModelsSt.getPixels(j);
        for (int i=0; i<lattice.length; ++i)
		{
			//int q = i % modelsize;
			//int r = i / modelsize;
			if (!hex_mask[i])//(!iswithinHex(q,r, modelsize/2))
			{	lattice[i] = 0;}
		}
        
        //ModelsSt = Models.getStack();
        while (ModelsSt.getSize() < modelnum)
		{
			//Models.setSlice( Models.getStackSize() );
			//IJ.run(Models,"Add Slice","");
			ModelsSt.addSlice("", new short[modelsize * modelsize]); //maybe modelarea
			if(dead[ModelsSt.getSize()-1] != ModelState.RESET)
			{
				IJ.log("WARNING: model" + (ModelsSt.getSize()) + " had a non Reset entry in weights.txt but was missing from " + Models.getTitle());
				dead[ModelsSt.getSize()-1] = ModelState.RESET; // a new frame has to bee initalized( reseted)
			}
		}
		Models.setStack(ModelsSt);
		
		DiffSt = Diffimg.getStack();
		while(DiffSt.getSize() < modelnum)
		{
			//Diffimg.setSlice( Diffimg.getStackSize() );
			//IJ.run(Models,"Add Slice","");
			DiffSt.addSlice("", new short[modelsize * modelsize]); //maybe modelarea
		}
        Diffimg.setStack(DiffSt);
        //DiffSt = Diffimg.getStack();
        //ModelsSt = Models.getStack();
        for(int m = 0; m < modelnum; ++m )
        {
			short[] pixels = (short[])ModelsSt.getPixels(m+1);
			short[] diffpx = (short[])DiffSt.getPixels(m+1);

			for (int i=0; i<diffpx.length; ++i)
			{
				diffpx[i] -= 32768;
			}
			int swap_count = 0;
			switch(dead[m])
			{
				case RESET: //"Reset"
				for (int i=0; i<pixels.length; ++i)
					{
						//int q = i % modelsize;
						//int r = i / modelsize;
						if (hex_mask[i])//(iswithinHex(q,r, modelsize/2))
						{
							pixels[i] = lattice[i];
						}
						diffpx[i] = (short)0; // aka latticediff[i]
					}
				
				case NOISIFY:
					if( (modelactive > 0) && noisy_init_models )
					{
						for (int i=0; i<pixels.length; ++i)
						{
							
							if (hex_mask[i])//(iswithinHex(q,r, modelsize/2))
							{
								int q = i % modelsize;
								int r = i / modelsize;
								if(	distToCenter(q,r) <= active_radius &&
									random.nextDouble() < noisification_level )
								{

									short pixold = pixels[i];
									int dx2 = -active_radius + random.nextInt(modelactive);
									int dzmin = -active_radius;
									int dzmax = active_radius - dx2;
									if (dx2 < 0)
									{
										dzmin = -active_radius - dx2;
										dzmax = active_radius;
									}
									int dz2 = dzmin + random.nextInt(dzmax-dzmin);

									int x2 = cx + dx2;
									int z2 = cz + dz2;

									int q2 = x2 + z2 / 2;
									int k = q2 + z2 * modelsize;

									short otherpix = pixels[k];
										//randomly redistribute intensities of the two pixels
									if( (pixold + otherpix) > 0)
									{
										//short swap = (short) random.nextInt(pixold + otherpix);
										short swap = (short) (pixold - 3 + random.nextInt(5)); 
										
										pixels[i] = swap;
										pixels[k] = (short) (pixold + otherpix - swap);

										diffpx[i] += (pixels[i] - pixold);
										diffpx[k] += (pixels[k] - otherpix);

										++swap_count;
									}
								}
							}

						}
					}
					dead[m] = ModelState.CONTINUE; //reset only once
				//rescale after adding noise
				case LATTICE:   //"Lattice"
				case CONTINUE: //"Continue"
				case FINISHED: //"Finished"
				case FIXED: //"Fixed"
				break;
			}
			//move back the diff values
			for (int i=0; i<diffpx.length; ++i)
			{
				diffpx[i] += 32768;
			}
		}
		
		double scaling = 1.0; 
		if(normalize_model_intensities)
		{
			normalize_model_intensities = false;
			double weighted_avg = 0.0;
			for(int m = 0; m < modelnum; ++m )
			{
				short[] pixels = (short[])ModelsSt.getPixels(m+1);
				int pixsum = 0;
				int pixcount = 0;
				for(int i=0; i<pixels.length; ++i)
				{
					//int q = i % modelsize;
					//int r = i / modelsize;
					if(hex_mask[i])//( iswithinHex(q,r, modelsize/2))
					{
						pixsum += pixels[i];
						++pixcount;
					}
				}
				weighted_avg += ((double)pixsum * modelweight[m])/(double)pixcount;
			}
			scaling = modelmean/weighted_avg;
		
			for(int m = 0; m < modelnum; ++m )
			{
				short[] pixels = (short[])ModelsSt.getPixels(m+1);
				short[] diffpx = (short[])DiffSt.getPixels(m+1);
				//signed unsigned HACK
				for (int i=0; i<diffpx.length; ++i)
				{
					diffpx[i] -= 32768;
				}
				for (int i=0; i<pixels.length; ++i)
				{
					int q = i % modelsize;
					int r = i / modelsize;
					if(hex_mask[i])//( iswithinHex(q,r, modelsize/2))
					{
						pixels[i] = (short)Math.round( scaling * pixels[i] );
						if(pixels[i] < 1) {pixels[i] = (short)1;}
						if(pixels[i] >= modelGraylevels) {pixels[i] = (short)(modelGraylevels - 1);}
						diffpx[i] = (short)Math.round( scaling * diffpx[i] );
						checksum += pixels[i];
					}
					else
					{
						diffpx[i] = (short)0;
						pixels[i] = (short)0;
					}
				}
				//undo signed unsigned HACK
				for (int i=0; i<diffpx.length; ++i)
				{
					diffpx[i] += 32768;
				}
			}
		
			IJ.log("Model scaling: " + format.format(scaling) + "   checksum: " + checksum);
		}
        
        while (ModelsSt.getSize() > modelnum)
        {
			ModelsSt.deleteLastSlice();
		}
		Models.setStack(ModelsSt);
		
		while (DiffSt.getSize() > modelnum)
        {
			DiffSt.deleteLastSlice();
		}
		Diffimg.setStack(DiffSt);
		Diff_Matcher matcher = null;
		if(sym_scaling)
		{
			matcher = new Diff_Matcher(Models, 1, modelnum, bondlength, (median_filtering?modelGraylevels:0) );
		}
		if(matcher == null)
		{
			Arrays.fill(sym_factors,1.0);
		}
		else
		{	
			for(int m = 0; m < modelnum; ++m)
			{
				sym_factors[m] = median_filtering ? matcher.summed_matches[m] : matcher.inv_avg_matches[m];	
			}	
		}
    }

	int advancePixel()
	{
		//System.out.println("advancing Pixel");
		if(blockvolume == 1)
		{
			if (incrs[0] == 3 ) {	return 1;}
			else return 0;//(force_wrap?1:0);
		}

		int skip_max = blockvolume;

		int[] bx = new int[1];
		int[] by = new int[1];
		int[] bm = new int[1];
		int lx = blocksize;
		int ly = blocksize;
		int lm = blockdepth;

		//original dynamic coordinates
		bx[0] = blockx;
		by[0] = blocky;
		bm[0] = blockm;

		//transformed dynamic coordinates
		int[] bx2 = bx;
		int[] by2 = by;
		int[] bm2 = bm;

		if(flip_block)
		{	//exchange transformed x and y
			int[] tmp = bx2;
			bx2 = by2;
			by2 = tmp;

			int tl = lx;
			lx = ly;
			ly = tl;
		}

		if(blockdepth_priority)
		{	//cylic shift of tranformed coordinates
			int[] tmp = bm2;
			bm2 = by2;
			by2 = bx2;
			bx2 = tmp;

			int tl = lm;
			lm = ly;
			ly = lx;
			lx = tl;
		}

		//transformed dynamic index
		int ind3d = (bx2[0] + by2[0] * lx  + bm2[0] * lx * ly ) ;
		int areaxy = lx * ly;
		do
		{
			--skip_max;
			if(reverse_block_depth)
			{	ind3d = (ind3d - 1)  % blockvolume;}
			else
			{	ind3d = (ind3d + 1)  % blockvolume;}
			if(ind3d < 0)
			{	ind3d += blockvolume;}

			bm2[0] = ind3d / areaxy;
			int ind2d = ind3d % areaxy;
			by2[0] = ind2d / lx;
			bx2[0] = ind2d % lx;
			//System.out.println("blockx,bx2,lx " + blockx + "," + bx[0] + "," + lx);
			//System.out.println("blocky,by2,ly " + blocky + "," + by[0] + "," + ly);
			//System.out.println("blockm,bm2,lm " + blockm + "," + bm[0] + "," + lm);
			//check against original dynamic index
		} while ( ( (incrs[bx[0] + by[0] * blocksize + bm[0] * blockarea] == 3) && (skip_max >= 0) )  );

		//get the original dynamic coordinates
		blockx = bx[0];
		blocky = by[0];
		blockm = bm[0];

		//System.out.println("to block x,y,m " + blockx + "," + blocky + "," + blockm + " :" + skip_max);

		if(skip_max >= 0)
		{   //we are still on the same block
			return 0; //(force_wrap?1:0);
		} 
		else
		{	//the block is finished
			return 1;
		} 

	}

	int advanceModel()
	{
		//System.out.println("advancing model");
		homem += blockdepth;
		if(homem >= modelnum)
		{
			homem -= modelnum;
			//System.out.println(" " + homem + " : " + skip_max);
			return 1;
		}
		if(shuffle_update)
		{
			int blockind = get_and_shuffle(blockid, permutation[homem]);
			//permutator.getIndex((blockid + homem * permutator.getRange()/modelnum) % permutator.getRange());
			int dim = modelsize/blocksize;
			homex = blocksize * (blockind % dim) + ( ((wrapnum % 2) == 1)?blocksize/2:0 );
			homey = blocksize * (blockind / dim) + ( ((wrapnum % 2) == 1)?blocksize/2:0 );
		}
		
		
		//System.out.println(" " + homem + " : " + skip_max);
		return 0; //(force_wrap?1:0); //the next model
	}

	void reportWrap()
	{
		//force_wrap = false;
		IJ.log("Wrap: " + (++wrapnum) + "/" + target_wraps + "   Pval= " + format.format(highestPval));
		++new_wraps;
		IJ.log("Finished Models: " + deadModelCount + "/" + target_models);
		IJ.log("Pixels active/total " + num_pixels_active + "/" + num_pixels_total);
		if(!rough_optimization )
		{
			IJ.log("bisects/quadfits : " + dbg_bisect + "/" + dbg_quadfit);
		}
		IJ.log("pixel tests: " + pixel_tests);
		IJ.log("moves good/total: " + dbg_good_moves + "/" + dbg_moves);
		long now = (System.currentTimeMillis());
		date.setTime(System.currentTimeMillis());
		IJ.log( dateFormat.format(date) + " Elapsed time: " + format.format((0.001*(now-timeoffset))) + "s");
		IJ.log("Tests per s: " + format.format( (double)(pixel_tests + dbg_weights)/(0.001*(now-timeoffset)) ) );
		
		if(pixel_tests > 0)
		{	
			writehist = true;
			//make only 1 step if switching or at cases
			tell_use_weights = (use_weights != next_use_weights); 
			use_weights = next_use_weights;
			tell_use_max = (use_max != next_use_max);
			use_max = next_use_max;
			weight_trials = tell_use_weights? 1 : (use_weights ? (4 * modelnum) : 1);
								//|| (wrapnum <= modelnum / multiActivations)  //works as intended but does not help with overshooting intensities at startup
			if( (const_weights  	) && !tell_use_weights) //no need for weight trials except for polling
			{	weight_trials = 0;}
			init_weights = true;
			request_cache_report = (debug_level > 0);
			request_autosave = (debug_level > 1);
			///This is always 0.0 at 17.09.2015
			/*
			if(sym_scaling)
			{ 
				matcher.update_matches();
				fresh_symmetries = true;
				System.out.println("double checking symmetry scaling factors");
				for(int v = 0; v < modelnum; ++v)
				{
					System.out.print("\t" + (matcher.inv_avg_matches[v]-sym_factors[v]) );
				}
				System.out.println();
			}
			*/ 
		}
		/*else
		{
			throw new RuntimeException(Macro.MACRO_CANCELED);
		}*/
		num_pixels_active = 0;
		num_pixels_total = 0;
		pixel_tests = 0;
		dbg_moves = 0;
		dbg_good_moves = 0;
		dbg_bisect = 0;
		dbg_quadfit = 0;
		//weights/cases will only be reported after adjusting them 
		
		if(alternating_scan_direction)
		{
			flip_model = !flip_model;
			flip_block = !flip_block;
		}
		
		timeoffset = now;
	}

	int advanceBlock()
    {
		//System.out.println("advancing block");
		if(flip_model) //exchange x and y before advancing
		{
			int tmp = homex;
			homex = homey;
			homey = tmp;
		}
		int reply = 0;
		boolean rerun = false;
		//if(force_wrap)
		//{
		//	blockid = -1;			
		//}
		do
		{
			rerun = false;
			if(spiral)
			{
				if(reply == 1)
				{
					homex = modelsize/2-blocksize + ( ((wrapnum % 2) == 1)?blocksize/2:0 );
					homey = modelsize/2-blocksize + ( ((wrapnum % 2) == 1)?blocksize/2:0 );
					spiralcountmax = 1;
					spiralDelta = 1;
					spiralcount = 0;
					spiraldirection = 0;
				}
				//go outward in a spiral
				if(spiralcount >= spiralcountmax)
				{
					//extend the range of the next horizontal row
					if( (spiraldirection == 1) || (spiraldirection == 3) )
					{
						spiralcountmax += spiralDelta;
					}
					spiraldirection = (spiraldirection + 1)%4; //turn clockwise
					spiralcount = 0; //reset counter
				}
				++spiralcount;
				switch(spiraldirection)
				{
					case 0: homex += blocksize; break;
					case 1: homey += blocksize; break;
					case 2: homex -= blocksize; break;
					case 3: homey -= blocksize; break;
					//default:
					//throw new IllegalStateException("spiral direction must be within [0,3]");
				}
				++blockid;
				if( (homex < 0) || (homex >= modelsize) ||
				(homey < 0) || (homey >= modelsize) || (spiralcountmax <= 0))
				{
					++reply;
					rerun = true;
					blockid = -1;
				}
			}
			else// shuffle_update || line scanning
			{
				int num_blocks = (modelsize/blocksize) * (modelsize/blocksize);
				if(++blockid >= num_blocks )
				{
					++reply;
					rerun = true;
					blockid = -1;
					//if(shuffle_update) //no longer needed since get_and_shuffle is used
					//{	shuffle(permutation);} //permutator.shuffle();
				}
				else
				{
					int blockind = shuffle_update? get_and_shuffle( blockid, permutation[homem]) : blockid; //permutator.getIndex(blockid)
					int dim = modelsize/blocksize;
					homex = blocksize * (blockind % dim) + ( ((wrapnum % 2) == 1)?blocksize/2:0 );
					homey = blocksize * (blockind / dim) + ( ((wrapnum % 2) == 1)?blocksize/2:0 );
				}	
			}
			/*if(force_wrap)
			{
				reply = 1;
				//blockid was already reinitialized
			}*/
			
			
			if(reply == 1)
			{
				//check if the scan pattern has changed
				spiral = scan_pattern.equals("spiral");
				shuffle_update = scan_pattern.equals("noise");
			}
		}
		while(rerun);
		
		if(reply > 1)
		{	IJ.log("ERROR: illegal reply in advanceBlock(): " + reply);}
		
		if(reply == 1)
		{
			if( wrapnum >= 0 ) // || force_wrap
			{
				reportWrap();
			}
			else //ignore and hide the first half wrap
			{
				++wrapnum;
				reply = 0; 
			}		
		}
		// !<= will also be true if highestPval_at_blockinit is NaN
		// but spare the effort if there was no increase in Pval anyways
		else if(use_weights && aggressive_weights && (!(highestPval <= highestPval_at_blockinit)) )
		{
			weight_trials = (use_weights ? (modelnum) : 1);
			init_weights = true;
			highestPval_at_blockinit = highestPval;
		}

		if(flip_model) //and change back after advancing
		{
			int tmp = homex;
			homex = homey;
			homey = tmp;
		}
		return reply;
	}

	void initWeightAdjustment(boolean trial)
	{
		/*
		for(int m=0; m<modelnum; ++m)
		{
			if(modelweight[m] != modelweightupd[m])
			System.out.println("Model" + m + " " + modelweight[m] + " != " + modelweightupd[m]);
		}
		*/
		int skip_max = modelnum;
		do
		{
			weightupdnum = (weightupdnum + 1) % (modelnum);
			if(--skip_max < 0) 
			{
				date.setTime(System.currentTimeMillis());
				System.out.println( dateFormat.format(date) );
				throw new IllegalStateException("no more active models");
			}
		} while ( dead[weightupdnum] == ModelState.FIXED ||  dead[weightupdnum] == ModelState.LATTICE); //skip the lattice and "Fixed" models

		oldwght = modelweight[weightupdnum];
		bestwght = oldwght;
		newwght = oldwght;
		if (trial)
		{
			newwght = oldwght * (0.9d+0.22d*random.nextDouble());
			if (newwght < 0.0001) {	newwght = 0.0001;}
			if (newwght > 0.9999) {	newwght = 0.9999;}
		}
		
		weight_incrs = -1;
		bestPvalW = highestPval;
		secondbestPvalW = highestPval;
		thirdbestPvalW = highestPval;
		//System.out.println("initWeighttrial\tm = " + weightupdnum + " w_incrs: " + weight_incrs + "\tbest: " + bestwght + "\tnew: " + newwght);
	}

	void weightAjustment()
	{
		if(weight_incrs >= 0)
		{
			double tmp = newwght;
			switch(weight_incrs)
			{
				case 2:
					if (reply > 0)
					{	newwght = (tmp + bestwght)/2;}
					else
					{	newwght = bestwght + 0.5*(bestwght - tmp);}
					break;
				case 1: //accelerated linear search
					if( reply > 0 )
					{
						newwght += 2*wstep;
					}
					else
					{
						newwght -= wstep/2;
						weight_incrs = 2;
					}
					break;
				case 0: //the first real step
					if(reply > 0)
					{
						newwght += 2*wstep;
						weight_incrs = 1;
					}
					else if (reply < 0)
					{
						newwght -= (3*wstep);
						weight_incrs = 1;
					}
					else //reply == 0
					{
						newwght -= wstep/2;
						weight_incrs = 2;
					}
					break;
			}

			oldwght = tmp;

		}
		else
		{ weight_incrs = 0;}
		if (newwght < 0.0001)
		{
			newwght = 0.0001;
		}
		if (newwght > 0.9999) {	newwght = 0.9999;}
		
		modelweightupd[weightupdnum ] = newwght;
		normalizeWeights(modelweightupd);
		newwght = modelweightupd[weightupdnum ];
		wstep = newwght-oldwght;
		//System.out.println("Weighttrial\tm = " + weightupdnum + " w_incrs: " + weight_incrs + "\tbest: " + bestwght + "\tnew: " + newwght);
		if(Math.abs(1-newwght/oldwght) < 0.005 )
		{
			//System.out.println("Weights converged m = " + weightupdnum + "\tdelta: " + Math.abs(1 - newwght/oldwght));
			renew = 1;
			if (Math.abs(1 - newwght/oldwght) < 0.001 )
			{	cmsok = 0; }
		}
		if(cmsok == 1) {++dbg_weights;}
	}

    int initadjustmentBlock(ImagePlus Models, ImagePlus Activity)
    {
		//System.out.println("initalizing Block");
		ongoing_Adj = 0;
		ImageStack ModelsSt = Models.getStack();
		ImageStack ActivitySt = Activity.getStack();
		//System.out.println("Initializing Block home x,y,m :" + homex + "," + homey + "," + homem);
		testing_avg = false;
		for(int bm = 0; bm < blockdepth; ++bm)
		{
			int mm = (homem+bm) % modelnum;
			short[] pixels = (short[])ModelsSt.getPixels( mm + 1 );
			short[] activity_pixels = (short[])ActivitySt.getPixels( mm + 1 );
			for(int ind = 0; ind < blockarea; ++ind)
			{
				int volind = ind + blockarea * bm;
				int by = (homey + ind / blocksize);
				int bx = (homex + ind % blocksize);
				// wrap onto model
				bx = ( modelsize + bx ) % modelsize;
				by = ( modelsize + by ) % modelsize;
				int bind = bx + modelsize * by;
				short val = (short)(activity_pixels[bind] - 32768);
				val /= 2; //exponential fading of activity
				ModelState ms = dead[ mm ];
				boolean inside_hexarea = hex_mask[bx+by*modelsize];//iswithinHex(bx, by, modelsize/2);
				if(inside_hexarea && (ms != ModelState.FINISHED) && (ms != ModelState.LATTICE) )
				{
					++num_pixels_total;
				}
				//FIXME the last test is actually need in rare cases, but why!?
				if ( (ms == ModelState.FINISHED) || (ms == ModelState.LATTICE) || 
					 (val < 0) || (!inside_hexarea) ) //|| ( (wrapnum < 1) && (mm == 0) )
				{	
					incrs[volind] = 3;//skip them
				} 
				else
				{
					++num_pixels_active;
					val -= (28 + random.nextInt(12)); //cost for picking pixel
					incrs[volind] = ( (val > 24) || (random.nextDouble()<0.05) )?-2:-1; //-1..normal tuning, -2 one extra move attempt
					++ongoing_Adj;
				}
				//initialize adjustment
				activity_pixels[bind] = (short)(val + 32768);
				oldvalinit[volind] = pixels[ bind ];
				bestval[volind] = oldvalinit[volind];
				oldval[volind] = oldvalinit[volind];
				secondbestval[volind] = oldvalinit[volind];
				thirdbestval[volind] = oldvalinit[volind];

				bestPval[volind] = highestPval; //the highest Pval
				secondbestPval[volind] =  highestPval; //the secondhighest Pval
				thirdbestPval[volind] =  highestPval; //the thirdhighest Pval
				int x = bx - (by-(by&1))/2;
				int z = by;
				int best_guess = 0;
				int num_guesses = 0;
				
				for(int dx = -1; dx <= 1; ++dx)
				{
					for(int dz = -1; dz <= 1; ++dz)
					{
						int dy = -dx - dz;
						if(Math.abs(dx)+Math.abs(dy)+Math.abs(dz) <= 2)
						{
							int nx = x + dx;
							int nz = z + dz;
							int ny = -nx -nz;
							int[] cube = {nx-cx,ny-cy,nz-cz};
							hexagonaltrap(cube, modelsize/2);
							nx = cube[0] + cx;
							ny = cube[1] + cy;
							nz = cube[2] + cz;
							int nq = nx + (nz-(nz&1))/2;
							int nr = nz;
							int nind = nq + nr * modelsize;
							best_guess += pixels[ nind ];
							if(		
								(incrs[volind] == -1) && 
								( 
									(wrapnum >= modelnum/multiActivations)  || 
									(distToCenter(nq, nr) < modelactive/2) 
								)
							  )
							{	activity_pixels[nind] -= 2;} //cost for picking a neighbor
							++num_guesses;
						}	
					}
				}
				best_guess = (int) Math.round((double)best_guess/(double)num_guesses);
				avg_blockval[volind] = best_guess;
				testing_avg |= (best_guess != oldvalinit[volind]);
				
				if(best_guess != oldvalinit[volind])
				{
					testval[volind] = best_guess; 
				}
				else
				{
					int bstep =	1 + random.nextInt(5);
					if(random.nextInt(2) == 1) {bstep = -bstep;}
					testval[volind] = best_guess + bstep;
					if(testval[volind] == oldvalinit[volind]) 
					{
						testval[volind] = oldvalinit[volind] + bstep;
					}
				}
				if(testval[volind] >= modelGraylevels - 1) { testval[volind] =  oldvalinit[volind] - 1;}
				if(testval[volind] <= 0) { testval[volind] = oldvalinit[volind] + 1;}
			}
		}
		//initialize coordinates at last pixel
		blockx = blocksize -1;
		blocky = blocksize -1;
		blockm = blockdepth -1;
		advancePixel();//move to first valid pixel
		if(ongoing_Adj > 0)
		{
			int num_blocks = (modelsize/blocksize) * (modelsize/blocksize);
			IJ.showProgress(blockid,num_blocks);
		}
		
		return 1;
	}


    int mktestadjustment3(int[] chx, int[] chy, int[] chm,
    int[] newv, int currentval, int chlistmax, int chlistlen)
    {
        
        int ind = blockx + blocky * blocksize + blockm * blockarea;
		if(incrs[ind] == -2)
		{
			incrs[ind] = -1;
			
			int[] oldv = new int[chx.length];
			chm[chlistlen] = (homem + blockm + modelnum) % modelnum;
			chx[chlistlen] = (homex + blockx + modelsize) % modelsize;
			chy[chlistlen] = (homey + blocky + modelsize) % modelsize;
			++chlistlen;
			
			int chlen = switch_pixels(Models, chx, chy, chm, newv, oldv);
			if(chlen == 2) //an actual switch can be attempted
			{	
				//the whole adjustment block will be dismissed if this attempt will be sucessfull
				strategy = Strategy.PIX_MOVES;
				++dbg_moves;
				return chlen;
			}
			//otherwise just proceed with fine tuning
		}
		
		
		
		if(incrs[ind] == -1) //fresh but maybe already blobbed
		{
			if(currentval == testval[ind])
			{
				if(currentval == 0)
				{	testval[ind] = 1;}
				else if (currentval == modelGraylevels - 1)
				{	testval[ind] = modelGraylevels - 2;}
				else if(random.nextInt(2) == 1) 
				{	++testval[ind];}
				else
				{	--testval[ind];}
			}
		}
		
		/*
		if(incrs[ind] == 3)
		{
			throw new IllegalStateException("invalid incrs == 3 in mktestadjustment3");
		}
		*/
		if( incrs[ind] != -1)
		{

			int bstep = testval[ind] - oldval[ind];
			if(incrs[ind]!=-100) //the EM value was only suggested not actually testet
			{	oldval[ind] = testval[ind];}
			switch (incrs[ind])
			{
				case 2: //we have identifyed the range of the maximum
					{
						/*	we assume that the max is within or close to the three best guesses
						 *  take the maximum of a quadratic function a*x²+b*x+c = y | a < 0
						 * 	datapoints dxij is difference in x, dyij in y
						 *  s is for difference between squares.
						 */
						boolean fallback = !enable_quad_fits;
						if(!fallback)
						{
							int dx12 = bestval[ind] - secondbestval[ind];
							int dx23 = secondbestval[ind] - thirdbestval[ind];
							double dy12 = bestPval[ind] - secondbestPval[ind];
							double dy23 = secondbestPval[ind] - thirdbestPval[ind];

							fallback = ! ((dx12 != 0) && (dx23 != 0) && (dy12 > 0.0) && (dy23 > 0.0));

							if( !fallback )
							{
								int dx12s = bestval[ind]*bestval[ind] - secondbestval[ind]*secondbestval[ind];
								int dx23s = secondbestval[ind]*secondbestval[ind] - thirdbestval[ind]*thirdbestval[ind];

								double dy12s = bestPval[ind]*bestPval[ind] - secondbestPval[ind]*secondbestPval[ind];
								double dy23s = secondbestPval[ind]*secondbestPval[ind] - thirdbestPval[ind]*thirdbestPval[ind];

								double a = ( dy12 * dx23  - dy23 * dx12  ) / ( dx12s * dx23  - dx23s * dx12  );
								double b = ( dy12 * dx23s - dy23 * dx12s ) / ( dx12  * dx23s - dx23  * dx12s );

								if(a != 0) { testval[ind] = (int) Math.round( -b / (2*a) );}

								fallback = (a >= 0);//There is no local max!
							}
						}
						if (fallback) //fallback to dampened oscillatory behavior
						{
							if(reply < 0 )
							{
								testval[ind] = bestval[ind] + (bestval[ind]-testval[ind])/2; // check the opposite side
							}
							else
							{
								testval[ind] = (bestval[ind] + secondbestval[ind])/2;
							}
							++dbg_bisect;
						}
						else
						{
							//we do only the first quadratic fit in quick mode
							++dbg_quadfit;
						}
						pixel_done = quick_optimization;

					}
				break;
				case 1:
					if ((reply <= 0)) //actually <= would be sufficient
					{
						incrs[ind] = 2; //now we have lower and upper bound
						bstep = -bstep/2;
						pixel_done = rough_optimization;
					}
					else //still monotone behavior keep accelerating with factor 2
					{
						bstep = 2*bstep;
					}
					testval[ind] += bstep;
				break;
				case 0: //we are in the beginning
					if ((reply<0)) //exchange test and old and update bstep
					{
						incrs[ind] = 1;
						testval[ind] = oldvalinit[ind] +  2 * (oldvalinit[ind]-testval[ind]);
						oldval[ind] = oldvalinit[ind];
						bstep = testval[ind] - oldval[ind];

					}
					else if ( (reply == 0) && !use_max && use_weights) //ok that should  happen really rarely
					{
						if( use_beamprofile || (debug_level > 1) ) ///DEBUG silence these for EM mode
						{
							System.out.println("Freak event dPval == 0.0 " + oldvalinit[ind] + "->" + testval[ind]);
							System.out.println("at x,y,m " + ((homex + blockx ) % modelsize) + "," + ((homey + blocky ) % modelsize) + "," + ((homem + blockm) % modelnum) + " block index: " + ind);
						}
						incrs[ind] = 2; //now we have lower and upper bound
						testval[ind] = ( testval[ind] + oldvalinit[ind] )/2;
					}
					else //ok just march this way
					{
						incrs[ind] = 1;
						bstep = 2*bstep;
						testval[ind] += bstep;
					}
				break;
				case -100: //em_maximized value
					//System.out.println("mktestadjustment3: found an EM prepared update: " + currentval + " -> "+ testval[ind]);
					if(currentval != testval[ind])
					{
						incrs[ind] = -200;
						renew = 0;
						break;
					}
					else
					{
						//trigger the unsc mechanism here
						
					}
				case -200: //em_maximized and tested or identical anyways
					renew = 1; //We are finished
					--ongoing_Adj;
					cmsok = 0; //choose a new strategy immediately
					incrs[ind] = 3; //pixel is optimized
					//System.out.println("mktestadjustment3: finished an EM tested value: " + currentval + " -> " + testval[ind]);
				break;
				default:
					date.setTime(System.currentTimeMillis());
					System.out.println( dateFormat.format(date) );
					throw new IllegalStateException("unrecogniced incrs in mktestajustment3");
			}

			if(testval[ind] >= modelGraylevels)
			{
				testval[ind] = oldval[ind] + (modelGraylevels-oldval[ind])/2;
			}
			else if(testval[ind] <= 0)
			{
				testval[ind] = oldval[ind]/2;
			}

			if( (	   (testval[ind] == currentval)
					|| (testval[ind] == oldvalinit[ind])
					|| (testval[ind] == oldval[ind])
					|| (testval[ind] == bestval[ind])
					|| (testval[ind] == secondbestval[ind]
					|| testval[ind] == thirdbestval[ind])     ) /*&& (incrs[ind] != -200)*/ )
			{
				if(incrs[ind] == -200)
				{
					
					date.setTime(System.currentTimeMillis());
					System.out.println( dateFormat.format(date) );
					System.out.println("oops killed test before commiting it");
					System.out.println("currentval: " + currentval);
					System.out.println("testval: " + testval[ind]);
					System.out.println("oldvalinit: " + oldvalinit[ind]);
					System.out.println("oldval: " + oldval[ind]);
					System.out.println("besval: " + bestval[ind]);
					System.out.println("secondbestval: " + secondbestval[ind]);
					System.out.println("thirdbestval: " + thirdbestval[ind]);
					
					throw new IllegalStateException("shit happened");
				}
				
				renew = 1; //We maybe finished
				--ongoing_Adj;
				cmsok = 0; //choose a new strategy immediately
				incrs[ind] = 3; //pixel is optimized
				/*System.out.println("Pixel block x,y,m " + homex + "," + homey + "," + homem + " converged " +
				oldvalinit[ind] + " -> " + bestval[ind] + " : remaining " + ongoing_Adj);*/
			}
		}
		else// if (incrs[ind] == -1)
		{
			incrs[ind] = 0;  //switch to binary search
		}
        if(renew == 0) //do (initial or further)step
        {
			cmsok=1;
			newv[chlistlen] = testval[ind];
			chm[chlistlen] = (homem + blockm + modelnum) % modelnum;
			chx[chlistlen] = (homex + blockx + modelsize) % modelsize;
			chy[chlistlen] = (homey + blocky + modelsize) % modelsize;
			++chlistlen;
			/*
			if(incrs[ind] == -200)
			{
				int l = chlistlen-1;
				System.out.println("prepared chlist for EM,  m: " + chm[l] + ", x: " + chx[l] + ", y: " + chy[l] + ", testing:" + testval[l]);
			}*/
			
        }
        return chlistlen;
    }

	void updateWeightAdjustment()
	{
		if(Pval >= bestPvalW)
		{
			thirdbestPvalW = secondbestPvalW;
			secondbestPvalW = bestPvalW;
			bestPvalW = Pval;

			thirdbestwght = secondbestwght;
			secondbestwght = bestwght;
			bestwght = newwght;
		}
		else if (Pval >= secondbestPvalW)
		{
			thirdbestPvalW = secondbestPvalW;
			secondbestPvalW = Pval;

			thirdbestwght = secondbestwght;
			secondbestwght = newwght;
		}
		else if (Pval >= thirdbestPvalW)
		{
			thirdbestPvalW = Pval;
			thirdbestwght = newwght;
		}
		else if (bestPvalW == secondbestPvalW) //set second best anyways
		{
			thirdbestPvalW = secondbestPvalW;
			secondbestPvalW = Pval;

			thirdbestwght = secondbestwght;
			secondbestwght = newwght;
		}
	}



	int updatetestajustment(ImagePlus Models)
	{
		int ind = blockx + blocky * blocksize + blockm * blockarea;
		double dPval = Pval - bestPval[ind];

		if (pixel_done)
		{
			pixel_done = false;
			renew = 1; //We maybe finished
			--ongoing_Adj;
			cmsok = 0; //choose a new strategy immediately //TODO maybe not needed here
			incrs[ind] = 3; //pixel is optimized
			return 0;
		}

		if(( dPval >= 0) ) //update our best guess so far for current ajustment
		{
			thirdbestval[ind] = secondbestval[ind];
			secondbestval[ind] = bestval[ind];
			bestval[ind] = testval[ind];
			thirdbestPval[ind] = secondbestPval[ind];
			secondbestPval[ind] = bestPval[ind];
			bestPval[ind] = Pval;
			
			ImageStack ModelsSt=Models.getStack();
			for(int sm = 0; sm < blockdepth; ++sm)
			{	
				short[] pixels = (short[])ModelsSt.getPixels((homem+sm+modelnum)%modelnum+1);
				for(int sy = 0; sy < blocksize; ++sy)
				{
					for(int sx = 0; sx < blocksize; ++sx)
					{
						int sind = sx + sy * blocksize + sm * blockarea;
						int mind = ( (modelsize + homex + sx) % modelsize ) + 
								   ( (modelsize + homey + sy) % modelsize ) * modelsize;
						if(sind != ind)//we are on a neighboring pixel that may have been "blobbed"
						{
							int dv = pixels[mind] - bestval[sind];
							if(incrs[sind] == -1) //not yet updated 
							{
								if(dv != 0) //we were actually "blobbed"
								{
									testval[sind] += dv;
									//Check to stay inside grayvalues
									if(testval[sind] >= modelGraylevels)
									{	
										testval[sind] = modelGraylevels -1;
										if(testval[sind] == pixels[mind])
										{
											--testval[sind];
										}	
									}
									else if(testval[sind] < 0)
									{
										testval[sind] = 0;
										if(testval[sind] == pixels[mind])
										{
											++testval[sind];
										}	
									}
									oldvalinit[sind] = pixels[mind];
									bestval[sind] = oldvalinit[sind];
									oldval[sind] = oldvalinit[sind];
									secondbestval[sind] = oldvalinit[sind];
									thirdbestval[sind] = oldvalinit[sind];
								}
								//update the starting Pval
								bestPval[sind] = Pval; //the highest Pval
								secondbestPval[sind] =  Pval; //the secondhighest Pval
								thirdbestPval[sind] = Pval; //the thirdhighest Pval
									
							}
							else if(incrs[sind] < 3) // already modified and still active
							{
								if(dv != 0) //avdvance as if it was a real update for that pixels
								{
									testval[sind] += dv; //shift the next sheduled test
									if(testval[sind] >= modelGraylevels)
									{	
										testval[sind] = modelGraylevels -1;
										if(testval[sind] == pixels[mind])
										{
											--testval[sind];
										}	
									}
									else if(testval[sind] < 0)
									{
										testval[sind] = 0;
										if(testval[sind] == pixels[mind])
										{
											++testval[sind];
										}	
									}
									
									
									
									thirdbestval[sind] = secondbestval[sind];
									secondbestval[sind] = bestval[sind];
									bestval[sind] += dv; //shift the latest value
									
									thirdbestPval[sind] = secondbestPval[sind];
									secondbestPval[sind] = bestPval[sind];
									bestPval[sind] += dPval; //the highest Pval
								} 
								else //simply translate pvals as if it would be uncorrelated
								{
									bestPval[sind] += dPval;
									secondbestPval[sind] += dPval;
									thirdbestPval[sind] += dPval;
								}
								
							} 
						}
					}
				}
			}
		}
		else if(Pval >= secondbestPval[ind]) //check for new second best
		{
			//System.out.println("new secondbest Pval " + testval[ind] + "\t" + Pval);
			thirdbestval[ind] = secondbestval[ind];
			secondbestval[ind] = testval[ind];
			thirdbestPval[ind] = secondbestPval[ind];
			secondbestPval[ind] = Pval;
		}
		else if( (Pval >= thirdbestPval[ind]) || (thirdbestPval[ind] == secondbestPval[ind]) ) //check for new third best Pval
		{
			//System.out.println("new thirdbest Pval " + testval[ind] + "\t" + Pval);
			thirdbestval[ind] = testval[ind];
			thirdbestPval[ind] = Pval;
		}
		else if (secondbestPval[ind] == bestPval[ind]) //initialize second anyways
		{
			thirdbestval[ind] = secondbestval[ind];
			secondbestval[ind] = testval[ind];
			thirdbestPval[ind] = secondbestPval[ind];
			secondbestPval[ind] = Pval;
		}
		return 0;
	}

	private boolean initialize_weights(String path, String Weighttitle)
	{
		boolean WeightstxtOk = true;
		if(!Weighttitle.equals("default") ) //attempt parsing file
		{
			//IJ.log("parsing Weights.txt");
			File weightfile = new File(path + Weighttitle);
			if(weightfile.isFile() && weightfile.setReadable(true) )
			{
				try
				{
					deadModelCount = 0;
					int j = 0;
					ModelState ded = ModelState.RESET;
					double weight = 0.0;
					Scanner sc = new Scanner(new FileReader(weightfile));
					//skip any header until a line starts with an value
					while(sc.hasNextLine() && !sc.hasNextInt() )
					{
						sc.nextLine();
					}
					for (int i=0; i < modelnum; ++i)
					{
						if(sc.hasNext() && sc.hasNextInt())
						{
							sc.nextInt(); //skip the index
							weight = sc.nextDouble();
							String state = sc.hasNextLine() ? sc.nextLine().trim() : modelStatus[0];
							System.out.println("" + weight + "\t" + state );
							if (state.equals(modelStatus[0]) || state.equals("") )
							{
								if(noisify_at_restart) {ded = ModelState.NOISIFY;}
								else {ded = ModelState.CONTINUE;}
								unsc[j] = 0;
							}
							else if (state.equals(modelStatus[1]) )
							{	ded = ModelState.FINISHED; ++deadModelCount;}
							else if (state.equals(modelStatus[2]) )
							{	ded = ModelState.FIXED; unsc[j] = 0;}
							else if (state.equals(modelStatus[3]) )
							{	ded = ModelState.LATTICE; ++deadModelCount;}
							else if (state.equals(modelStatus[4]) )
							{	ded = ModelState.RESET; unsc[j] = 0;}
							else if (state.equals(modelStatus[5]) )
							{	ded = ModelState.NOISIFY; unsc[j] = 0;}
							else
							{
								System.out.println("Unrecognized modelstate\"" + state + "\" Please fix Weights.txt. ");
								IJ.log("Unrecogniced modelstate +\"" + state + "\"+ in Weight.txt");
								WeightstxtOk = false;
								return false;
							}
						}
						else
						{
							ded = ModelState.RESET; //shall be reseted
							weight = default_weight;
							unsc[j] = 0;
						}
						if(weight >= 0.0) //might actually be zero if cases are used
						{
							modelweight[j] = weight;
							dead[j] = ded;
							++j;
						}
					}
					sc.close();
				}
				catch(Exception e)
				{
					e.printStackTrace();
					WeightstxtOk = false;
				}
			}
			else
			{
				WeightstxtOk = false;
			}

		}

		if(Weighttitle.equals("default") || !WeightstxtOk) //saftey fallback to default weights
		{
			deadModelCount = 0;
			if(!WeightstxtOk)
			{
				IJ.log("missing or corrupted \"Weights.txt\", assuming default weights and model states");
				//IJ.log("default model states will reset all but the first model");
			}

			modelweight[0] = (modelnum>1)? init_lattice_weight : 1.0;
			/*if(do_full_EM)
			{
				dead[0] = ModelState.NOISIFY;
			}
			else*/
			{
				dead[0] = (do_pix_EM || do_full_EM) ? ModelState.RESET : ModelState.LATTICE;
				++deadModelCount;
			}
			//System.out.println("" + (0+1) + "\t" + modelweight[0] + "\t" + getModelStateString(dead[0]));
			double sofar = modelweight[0];
			for (int i=1; i<modelnum; ++i)
			{
				//classic equal weight distribution
				if(fair_weights)
				{
					modelweight[i] = (1.0d-init_lattice_weight)/((double)(modelnum-1));
				}
				else
				{
					//unfair weight distribution
					if(i < modelnum-1 )
					{
						
						double fair_share = (1.0-sofar)/(modelnum-i);
						//weights should be tentatively growing so model 2 becomes
						//most flexible and heavier later models can compensate intensities
						//modelweight[i] = 0.25 * fair_share + 0.75 * random.nextDouble()*( fair_share );
						
						
						//Weights are tentatively shrinking so the last model becomes the most flexible
						//it is also activated as first model
						modelweight[i] = fair_share * ( 1.0 + 0.8 * random.nextDouble());
						
					}
					else
					{
						modelweight[i] = (1.0 - sofar);
					}
					sofar += modelweight[i];
				}
				dead[i] = ModelState.RESET;
			    //System.out.println("" + (i+1) + "\t" + modelweight[i] + "\t" + getModelStateString(dead[i]));
			}
		}

		normalizeWeights(modelweight); //in case the lattice weight has been altered manually
		//double weightSum = 0.0;
		for (int i=0; i<modelnum; ++i)
		{
			//weightSum += modelweight[i];
			modelweightupd[i] = modelweight[i];
		}
		//IJ.log("Sum of Weights: " + weightSum);
		return true;
	}

    // returns either a match or the last trial
    File locateNodes(String[] paths)
    {
		int i = 0;
		File nodes = new File (paths[i] + "nodes.txt");
		while(!(nodes.isFile() && nodes.setReadable(true)) && i < paths.length)
		{
			++i;
			nodes = new File (paths[i] + "nodes.txt");
		}
		return nodes;
	}

	//http://stackoverflow.com/questions/237159/whats-the-best-way-to-check
	//-to-see-if-a-string-represents-an-integer-in-java/237204#237204
	public boolean isInteger(String str) 
	{
		if (str == null) 
		{	return false;}
		int length = str.length();
		if (length == 0)
		{	return false;}
		int i = 0;
		if (str.charAt(0) == '-')
		{
			if (length == 1) 
			{	return false;}
			i = 1;
		}
		for (; i < length; i++)
		{
			char c = str.charAt(i);
			if (c <= '/' || c >= ':') 
			{	return false;}
		}
		return true;
	}

    int loadNodes(String path0, boolean reset_threads)
    {
		String[] places = {path0, IJ.getDirectory("imagej")};
		File nodes = locateNodes(places);
		if(nodes.isFile() && nodes.setReadable(true))
		{
			try
			{
				IJ.log(nodes.getPath());
				HashMap<String,Integer> busy_bindings = null;
				boolean has_busy = false;
				
				if(check_free_nodes)
				{
					IJ.log("scanning network for available workers ...");
					File busy_script = new File(IJ.getDirectory("imagej") + "bn.sh");
					if (busy_script.isFile() && busy_script.canExecute())
					{
						has_busy = true;
						busy_bindings = new HashMap<String,Integer>();
						java.lang.Runtime rt = java.lang.Runtime.getRuntime();
						Process ps = rt.exec( busy_script.getPath() );
						java.io.InputStream instr = ps.getInputStream();
						BufferedReader readstream = new java.io.BufferedReader(new InputStreamReader(instr));
						ps.waitFor(); //flood the Buffer with the script output
						String line = null;
						while( readstream.ready() )
						{
							line = readstream.readLine();
							Integer running = null;
							String host = line;
							if( readstream.ready() )
							{
								line = readstream.readLine();
								while(!isInteger(line) && readstream.ready())
								{	line = readstream.readLine();}
								
								running =  isInteger(line) ? new Integer( line ) : null;
							}
							//there may be several entries for one host but already the first reports on all running hexclients
							if ( !host.equals("<unknown>") && (running != null) && (running.intValue() > 0) && !busy_bindings.containsKey(host) )
							{
								//System.out.println("detected " + running + " workers at " + host);
								busy_bindings.put(host, running);
							}
						}
						ps.destroy(); //also closes all readers
					}
					else
					{
						IJ.log("missing shell script " + busy_script.getPath() );
					}					
				}
				else
				{
					IJ.log("Please make sure that these workers are available");
				}
				
				if(reset_threads)
				{
					threadVec.clear();
					threadadd = 0;
					use_threads = 0;
				}
				//we might only find less nodes, but that should not be an issue
				int j = 0;
				double weight = 0.0;
				int skip_nodes = check_free_nodes?0:threadadd; //threadadd; // zero should be fine if we scanned for already launched workers
				Scanner sc = new Scanner(new FileReader(nodes));
				//skip any header until a line starts with an value
				while(sc.hasNextLine() && !sc.hasNextInt() )
				{	sc.nextLine();}
				int i = 0;
				int next_rank = 0;
					
				while(sc.hasNext() && sc.hasNextInt())
				{
					int multiplicity = sc.nextInt();
					double perf = sc.nextDouble();
					double load = sc.nextDouble();
					boolean local = sc.nextBoolean();
					String executable = sc.nextLine().trim(); //gets whatever comes before newline
					String host = null;
					int port = 22;
					if(!local)
					{
						port = -1;
						Scanner sc2 = new Scanner(executable);
						sc2.next(); //ssh
						sc2.next(); //-p
						port = sc2.nextInt();
						host = sc2.next();
						sc2.close();
						
					}
					String server = local? "localhost" : (host + "_p" + port);					
					
					if( has_busy && busy_bindings.containsKey(server) && (port >= -1) )
					{
						Integer runInt = busy_bindings.get(server);
						int running = runInt.intValue();
						int more_running = running - multiplicity;
						if( more_running < 0 )
						{	more_running = 0;} 
						busy_bindings.put(server, new Integer(more_running));
						multiplicity -= running;
						//System.out.println("locking " + running + " workers at " + server);
						if(multiplicity < 0)
						{	multiplicity = 0;}
					}
					
					while( (multiplicity > 0) && (skip_nodes > 0))
					{
						--multiplicity;
						--skip_nodes;
					}
					if(multiplicity > 0)
					{
						addnode( executable, perf, load, local, next_rank, multiplicity);
						next_rank += multiplicity;
					}
				}	
				return next_rank;	
			}
			catch(Exception e)
			{
				e.printStackTrace();
				return -1;
			}	
		}
		else
		{
			IJ.log("could not open " + nodes.getPath() );
			return -1;
		}
	}
   
    void addnode(String next_command, double perf, double load, boolean local, int next_rank, int multiplicity)
    {
        int rank_end = next_rank + multiplicity;
        for (int i = next_rank; i < rank_end; ++i)
        {
            threadVec.add( new ThreadData() );
            threadVec.get(threadadd).command = next_command;
            threadVec.get(threadadd).can_write = local;
            threadVec.get(threadadd).performance = perf * load; //always use the product, except when writing a new nodes file
            threadVec.get(threadadd).load = load;
            threadVec.get(threadadd).worker = threadadd;
            ++threadadd;
        }
        IJ.log("" + multiplicity  + " workers: " + next_command + " performance: " + perf); 
    }
    
    boolean employThreads(int num_threads) //ideally called once per job
    {
		if(use_threads + num_threads > threadadd)
		{
			num_threads = threadadd - use_threads;
			IJ.log("Error: There are only " + num_threads + " workers available");
			return false;
		}
		
		if(num_threads < 1)
		{	return false;}
		
		for(int i = rank_id; i < rank_id + num_threads; ++i)
		{
			threadVec.get(use_threads).job = job_id;
			threadVec.get(use_threads).rank = i;
			++use_threads;
		}
		rank_id += num_threads;
		jobVec.get(job_id).performances = new double[rank_id];
		for (int i = 0; i < use_threads; ++i) //collect the performances
		{
			if(threadVec.get(i).job == job_id)
			{
				jobVec.get(job_id).performances[threadVec.get(i).rank] = threadVec.get(i).performance;
			}
		}
		return true;
	}
    
    boolean skipThreads(int num_threads) //ideally called once per job
    {
		if(use_threads + num_threads > threadadd)
		{
			num_threads = threadadd - use_threads;
			IJ.log("Error: There are only " + num_threads + " workers available");
			return false;
		}
		
		if(num_threads < 1)
		{	return false;}
		
		for(int i = 0; i < num_threads; ++i)
		{
			threadVec.remove(use_threads);
		}
		
		threadadd -= num_threads;
		max_threads -= num_threads;
		return true;
	}
    
    
    
	String getModelStateString(ModelState ms)
	{
		switch(ms)
		{
			case CONTINUE: return modelStatus[0];
			case FINISHED: return modelStatus[1];
			case FIXED:    return modelStatus[2];
			case LATTICE:  return modelStatus[3];
			case RESET:    return modelStatus[4];
			case NOISIFY:  return modelStatus[5];
			default:       return modelStatus[6];
		}
	}

	char getStrategyTag()
	{
		switch (strategy)
		{
			case PIXELADJ   :  return 'P';
			case WEIGHTADJ  :  return 'W';
			case PIXEL_EM   :  return 'p';
			case MODEL_EM	:  return 'm';	
			case TOTAL_EM	:  return 't';
			case POLL_HISTO :  return 'h';
			case PIX_MOVES  :  return 'M';
			case NOSTRATEGY :  return '?';
			default         :  return '*';
		}
	}


	private boolean loadOptions(File options, boolean master)
	{
		try
		{
			if(options.isFile() && options.setReadable(true) )
			{
				Scanner fs = new Scanner(new FileReader(options));
				while( fs.hasNextLine() )
				{
					String line = fs.nextLine();
					if(!line.equals(""))
					{
						Scanner ls = new Scanner(line);
						String token = ls.next();
						token = token.trim();
						if (token.startsWith("#"))
						{
							if(token.equals("#appname"))
							{	
								if(!getClass().getSimpleName().equals( ls.nextLine().trim() ))
								{
									IJ.error("invalid file " + options.getPath());
									return false;
								}
							}
							else if(token.equals("#sourcepath"))
							{	sourcepath = ls.nextLine().trim();}
							else if(token.equals("#Input_Data"))
							{	oldDatatitle = prefix + ls.nextLine().trim();}
							else if(token.equals("#Beam_Profile"))
							{	oldbeamProfileTitle = prefix + ls.nextLine().trim();}
							else if(token.equals("#noise_level"))
							{	noise_level = ls.nextDouble();}
							else if(token.equals("#Action"))
							{	if(master) Modus = ls.nextLine().trim();}
							else if(token.equals("#Models"))
							{	if(master) oldModeltitle = ls.nextLine().trim();}
							else if(token.equals("#Ptable"))
							{	oldPtableTitle = prefix + ls.nextLine().trim();}
							else if(token.equals("#Weights"))
							{}//{	if(master) Weighttitle = ls.nextLine().trim();}
							else if(token.equals("#modelnum"))
							{	if(master) modelnum = ls.nextInt();}
							else if(token.equals("#shadownum"))
							{	if(master) shadownum = ls.nextInt();}
							else if(token.equals("#parallel_activations"))
							{	if(master) multiActivations0 = ls.nextInt();}
							else if(token.equals("#modelactive"))
							{	if(master) modelactiveinit = ls.nextInt();}
							else if(token.equals("#bondlength"))
							{	if(master) bondlength = ls.nextInt();}
							else if(token.equals("#wobble"))
							{	if(master) wobble = ls.nextInt();}
							else if(token.equals("#modelGraylevels"))
							{	if(master) modelGraylevels = ls.nextInt();}
							else if(token.equals("#modelBrightness"))
							{	if(master) modelBrightness = ls.nextDouble();}
							else if(token.equals("#blocksize"))
							{	if(master) initblocksize = ls.nextInt();}
							else if(token.equals("#blockdepth"))
							{	if(master) initblockdepth = ls.nextInt();}
							else if(token.equals("#pvalscaling"))
							{	if(master) pvalscaling = ls.nextDouble();}
							else if(token.equals("#skipthreads"))
							{	if(master) skip_threads = ls.nextInt();}
							else if(token.equals("#threadnum"))
							{	initthreadnum = ls.nextInt();}
							else if(token.equals("#thread_offset"))
							{	 if(master) thread_offset = ls.nextInt();}
							else if(token.equals("#weightfr"))
							{	/*if(master) weightfr = ls.nextDouble();*/}
							else if(token.equals("#target_models"))
							{	if(master) init_target_models = ls.nextInt();}
							else if(token.equals("#target_wraps"))
							{	if(master) init_target_wraps = ls.nextInt();}
							else if(token.equals("#use_frames"))
							{	init_use_frames = ls.nextInt();}
							else if(token.equals("#pvaloffset"))
							{	pvaloffset = ls.nextDouble();}
							else if(token.equals("#external_merit"))
							{	external_merit = ls.nextBoolean();}
							else if(token.equals("#autoptableoffset"))
							{	autoptableoffset = ls.nextBoolean();}
							else if(token.equals("#initbenchmark"))
							{	if(master) initbenchmark = ls.nextBoolean();}
							else if(token.equals("#shuffle_frames"))
							{	shuffle_frames = ls.nextBoolean();}
							else if(token.equals("#spiral"))
							{	/*if(master) spiral = ls.nextBoolean();*/}
							else if(token.equals("#shuffle_update"))
							{	/*if(master) shuffle_update = ls.nextBoolean();*/}
							else if(token.equals("#scan_pattern"))
							{	if(master) scan_pattern = ls.nextLine().trim();}
							else if(token.equals("#alternating_scan_direction"))
							{	if(master) alternating_scan_direction = ls.nextBoolean();}
							else if(token.equals("#blockdepth_priority"))
							{	if(master) blockdepth_priority = ls.nextBoolean();}
							else if(token.equals("#correlated_optimization"))
							{	if(master) correlated_optimization = ls.nextBoolean();}
							else if(token.equals("#quick_optimization"))
							{	if(master) quick_optimization = ls.nextBoolean();}
							else if(token.equals("#rough_optimization"))
							{	if(master) rough_optimization = ls.nextBoolean();}
							else if(token.equals("#enable_quad_fits"))
							{	if(master) enable_quad_fits = ls.nextBoolean();}
							else if(token.equals("#noisify_at_restart"))
							{	if(master) noisify_at_restart = ls.nextBoolean();}
							else if(token.equals("#fixed_lattice_weight"))
							{	if(master) fixed_lattice_weight = ls.nextBoolean();}
							else if(token.equals("#frozen_lattice"))
							{	if(master) frozen_lattice = ls.nextBoolean();}
							else if(token.equals("#previous_seed"))
							{	if(master) previous_seed = ls.nextBoolean();}
							else if(token.equals("#sym_scaling"))
							{	if(master) sym_scaling = ls.nextBoolean();}
							else if(token.equals("#do_full_EM"))
							{	if(master) do_full_EM = ls.nextBoolean();}
							else if(token.equals("#do_pix_EM"))
							{	if(master) do_pix_EM = ls.nextBoolean();}
							else if(token.equals("#force_pix_EM"))
							{	if(master) force_pix_EM = ls.nextBoolean();}
							else if(token.equals("#median_filtering"))
							{	if(master) median_filtering = ls.nextBoolean();}
							else if(token.equals("#clone_models"))
							{	if(master) clone_models = ls.nextBoolean();}
							else if(token.equals("#reset_activity"))
							{	if(master) reset_activity = ls.nextBoolean();}
							else if(token.equals("#skip_first_clone"))
							{	if(master) skip_first_clone = ls.nextBoolean();}
							else if(token.equals("#interactive_cloning_session"))
							{	if(master) interactive_cloning_session = ls.nextBoolean();}
							else if(token.equals("#clone_weight"))
							{	if(master) clone_weight = ls.nextDouble();}
							else if(token.equals("#master_shows_dialog"))
							{	if(master) master_shows_dialog = ls.nextBoolean();}
							else if(token.equals("#collect_Histos"))
							{	if(master) collect_Histos = ls.nextBoolean();}
							else if(token.equals("#lattice_hopping"))
							{	if(master) lattice_hopping = ls.nextBoolean();}
							else if(token.equals("#looptotal"))
							{	if(master) looptotal = ls.nextInt();}
							else if(token.equals("#noisify_avg_data"))
							{	noisy_init_models = ls.nextBoolean();}
							else if(token.equals("#area_noise_level"))
							{	if(master) area_noise_level = ls.nextDouble();}
							else if(token.equals("#noisification_level"))
							{	if(master) noisification_level = ls.nextDouble();}
							else if(token.equals("#default_weight"))
							{	if(master) default_weight = ls.nextDouble();}
							else if(token.equals("#init_lattice_weight"))
							{	if(master) init_lattice_weight = ls.nextDouble();}
							else if(token.equals("#normalize_model_intensities"))
							{	if(master) normalize_model_intensities = ls.nextBoolean();}
							else if(token.equals("#fair_weights"))
							{	if(master) fair_weights = ls.nextBoolean();}
							else if(token.equals("#use_weights"))
							{	if(master) use_weights = ls.nextBoolean();}
							else if(token.equals("#const_weights"))
							{	if(master) const_weights = ls.nextBoolean();}
							else if(token.equals("#use_max"))
							{	if(master) use_max = ls.nextBoolean();}
							else if(token.equals("#black_init"))
							{	if(master) black_init = ls.nextBoolean();}
							else if(token.equals("#exclusive_autosave"))
							{	if(master) exclusive_autosave = ls.nextBoolean();}
							else if(token.equals("#keep_Idata"))
							{	if(master) keep_Idata = ls.nextBoolean();}
							else if(token.equals("#duplicate_Idata"))
							{	duplicate_Idata = ls.nextBoolean();}
							else if(token.equals("#real_log2"))
							{	if(master) real_log2 = ls.nextBoolean();}
							else if(token.equals("#check_free_nodes"))
							{	check_free_nodes = ls.nextBoolean();}
							else if(token.equals("#blobs"))
							{	if(master) blobs = ls.nextBoolean();} 
							else if(token.equals("#frame_shuffle_seed"))
							{	seed = ls.nextLong();}
							else if(token.equals("#report_interval"))
							{	if(master) report_interval = ls.nextInt();}
							else if(token.equals("#report_num"))
							{	if(master) report_num = ls.nextInt();}
							else if(token.equals("#report_cols"))
							{	if(master) report_cols = ls.nextInt();}
							else if(token.equals("#weight_cols"))
							{	if(master) weight_cols = ls.nextInt();}
							else if(token.equals("#debug_mode")) //backward compability conversion
							{	if(master) debug_level = (ls.nextBoolean()? 1 : 0) ;}
							else if(token.equals("#debug_level"))
							{	if(master) debug_level = ls.nextInt();}
							else if(token.equals("#another_job") || token.equals("#another_task")) //backward compability
							{	another_job = ls.nextBoolean();}
							else if(token.equals("#max_supported_workers"))
							{	if(master) MAX_THREADS = ls.nextInt();}
							else if(token.equals("#aggressive_weights"))
							{	aggressive_weights = ls.nextBoolean();}
							else if(master)
							{
								System.out.println("WARNING unknown token in options.txt: " + token + "\t" + (ls.hasNext()?ls.next():""));
							}

						}
						ls.close();
					}
				}
				fs.close();
				next_use_weights = use_weights;
				next_use_max = use_max;
				initpvaloffset = pvaloffset;
			}
			return true;
		}
		catch (Exception e)
		{
			System.out.println("Error: could not read " + options);
			e.printStackTrace();
			return false;
		}
	}

	void saveOptions(File options, boolean tipps)
	{
		try
		{
			java.io.FileOutputStream opts = new java.io.FileOutputStream(options, false);
			java.io.PrintWriter optp = new java.io.PrintWriter(opts);
			boolean master = (job_id == 0);
			optp.println("options are read from lines that start with #");
			optp.println("Please make sure your filenames do not contain any blank characters");
			if(tipps) optp.println("if the next option file should be included options.txt, options0.txt, options1.txt ...");
			optp.println("#appname\t" + getClass().getSimpleName());
			optp.println("#another_job\t" + another_job);
			if(tipps) optp.println("If images are not yet open they will be read from disk");
			if(tipps) optp.println("Input_Data will only be loaded from sourcepath");
			if(tipps) optp.println("Models and weights from source if they are not found in this directory");
			optp.println("#sourcepath\t" + sourcepath);
			int prfx_len = master_prefix.length();
			optp.println("#Input_Data\t" + Datatitle.substring(prfx_len));
			optp.println("#Beam_Profile\t" + beamProfileTitle.substring(prfx_len));
			optp.println("#Ptable\t" + master_PtableTitle.substring(prfx_len));
			if(tipps) optp.println("if greater zero Input_Data will be noisified so that the mean is divided by noise_level");
			optp.println("#noise_level\t" + noise_level);
			if(master) 
			{
				if(tipps) optp.println("external_merit provides pval as merit for e.g. Smart_Atomzer");
				optp.println("#external_merit\t" + external_merit);
				optp.println("#Action\t" + Modus);
				optp.println("#Models\t" + Modeltitle);
				//optp.println("#smooth_passes\t" + smooth_passes);
				//optp.println("#Weights\t" + Weighttitle);
				if(tipps) optp.println("The number of models may be changed at relaunching an optimization");
				optp.println("#modelnum\t" + modelnum);
				if(tipps) optp.println("Shadows are finely rotated within pi/3");
				optp.println("#shadownum\t" + shadownum);
				optp.println("#parallel_activations\t" + multiActivations0);
				if(tipps) optp.println("The diameter of the noisy hexagon in the center of a new model");
				optp.println("#modelactive\t" + modelactiveinit);
				if(tipps) optp.println("The bondlength in hexagonal pixels");
				optp.println("#bondlength\t" + bondlength);
				if(tipps) optp.println("wobble radius defines how close to lattice points translations are considered");
				optp.println("#wobble\t" + wobble);
				if(tipps) optp.println("Number of gray levels for models");
				optp.println("#modelGraylevels\t" + modelGraylevels);
				if(tipps) optp.println("Actual pixelvalue for total \"White\"");
				optp.println("#modelBrightness\t" + modelBrightness);
				if(tipps) optp.println("The pixeloptimization advances in entire blocks");
				if(tipps) optp.println("in every model the block is a square with edge length");
				optp.println("#blocksize\t" + initblocksize);
				if(tipps) optp.println("The depth of a block goes throug consecutive models");
				optp.println("#blockdepth\t" + initblockdepth);
				if(tipps) optp.println("pvalscaling is an effective power for summation over orbits");
				optp.println("#pvalscaling\t" + pvalscaling);
				if(tipps) optp.println("Number of workers to skip");
				optp.println("#skipthreads\t" + skip_threads);
				if(tipps) optp.println("lattice_hopping eenables or disables translational symmetries");
				optp.println("#lattice_hopping\t" + lattice_hopping);
				
			}
			if(tipps) optp.println("Number of workers to launch for reconstruction");
			optp.println("#threadnum\t" + initthreadnum);
			
			
			if(master)
			{ 
				if(tipps) optp.println("Offset to ID of workers, affects cpu binding");
				optp.println("#thread_offset\t" + thread_offset);
				//if(tipps) optp.println("Probability to do a weight optimization, the probability of optimizing");
				//if(tipps) optp.println("the next COMPLETE block is 1 - weightfr");
				//optp.println("#weightfr\t" + weightfr);
				if(tipps) optp.println("reconstruction finishes after reaching either target");
				optp.println("#target_models\t" + init_target_models);
				optp.println("#target_wraps\t" + init_target_wraps);
			}
			if(tipps) optp.println("Number of frames to use for reconstruction");
			optp.println("#use_frames\t" + init_use_frames);
			if(tipps) optp.println("Quit after first updated Pval");
			if(master) optp.println("#initbenchmark\t" + initbenchmark);
			if(tipps) optp.println("Only the first \"frames_used\" are shuffeled");
			optp.println("#shuffle_frames\t" + shuffle_frames);
			if(master)
			{
				//if(tipps) optp.println("If the new blocks are visited in a spiral pattern");
				//optp.println("#spiral\t" + spiral);
				//if(tipps) optp.println("The blocks will be visited in quasi random order");
				//optp.println("#shuffle_update\t" + shuffle_update);
				if(tipps) optp.println("Pattern for visiting blocks");
				optp.println("#scan_pattern\t" + scan_pattern);
				if(tipps) optp.println("normalize_model_intensities rescales the models");
				optp.println("#normalize_model_intensities\t" + normalize_model_intensities);
				if(tipps) optp.println("if the scann orientation (x <-> y) is switched every wrap");
				optp.println ("#alternating_scan_direction\t" + alternating_scan_direction);
				if(tipps) optp.println("noisify_avg_data adds some random noise to avg. data");
				optp.println("#noisify_avg_data\t" + noisy_init_models);
				if(tipps) optp.println("blockdepth is scanned before x and y");
				optp.println("#blockdepth_priority\t" + blockdepth_priority);
				if(tipps) optp.println("All pixels in a block are changed in turns");
				optp.println("#correlated_optimization\t" + correlated_optimization);
				if(tipps) optp.println("weights are scaled to compensate lowered model symmetries");
				optp.println("#sym_scaling\t" + sym_scaling);
				if(tipps) optp.println("employ global Expectation Maximation instead of PixelUpdates");
				optp.println("#do_full_EM\t" + do_full_EM);
				if(tipps) optp.println("sporadically attempt and check sinle pixel EM");
				optp.println("#do_pix_EM\t" + do_pix_EM);
				if(tipps) optp.println("allways do pixel EM and always accept the changes");
				optp.println("#force_pix_EM\t" + force_pix_EM);
				//if(tipps) optp.println("raise symmetry factors to this power");
				//optp.println("sym_fac_exponent\t" + sym_fac_exponent);
				if(tipps) optp.println("filter out lattice when calculating model symmetries (BAD IDEA!)");
				optp.println("#median_filtering\t" + median_filtering);
				if(tipps) optp.println("Debug tool attemping to collect the Ptable from the state of the reconstuction");
				optp.println("#collect_Histos\t" + collect_Histos);
				if(tipps) optp.println("Pixeladjustment stops after the second bisection or first quad fit");
				optp.println("#quick_optimization\t" + quick_optimization);
				if(tipps) optp.println("The \"lattice\" aka model0 is always the template for spawning new models");
				optp.println("#fixed_lattice_weight\t" + fixed_lattice_weight);
				if(tipps) optp.println("If model0 should be modified at all");
				optp.println("#frozen_lattice\t" + frozen_lattice);
				if(tipps) optp.println("the initial weight for model0 (regardless of its actual state)");
				optp.println("#init_lattice_weight\t" + init_lattice_weight);
				if(tipps) optp.println("weights are good for finding molecules, no weights are better for refinements");
				optp.println("#use_weights\t" + use_weights);
				if(tipps) optp.println("constant and equal weights do not discriminate the models");
				optp.println("#const_weights\t" + const_weights);
				if(tipps) optp.println("use_max (only best fit no averaging) is really no good for finding molecules, but superior for refinement");
				optp.println("#use_max\t" + use_max);
				if(tipps) optp.println("black initialization is much faster, but may run out of double range");
				optp.println("#black_init\t" + black_init);
				if(tipps) optp.println("with exclusive autosave and debuglevel > 1, no images are saved on halt, error or exit");
				optp.println("#exclusive_autosave\t" + exclusive_autosave);
				if(tipps) optp.println("keep Idata after recontruction has finished or jobfiles were written");
				optp.println("#keep_Idata\t" + keep_Idata);
				if(tipps) optp.println("worker will return a pval as if pvaloffset would have been zero");
				optp.println("#real_log2\t" + real_log2);
				if(tipps) optp.println("if models should be duplicated at a clean shut down");
				optp.println("#clone_models\t" + clone_models);
				if(tipps) optp.println("fraction of the weight beeing transfereed to the clones");
				optp.println("#clone_weight\t" + clone_weight);
				if(tipps) optp.println("if the first model (=lattice) should be skipped at cloning");
				optp.println("#skip_first_clone\t" + skip_first_clone);
				if(tipps) optp.println("interactive cloning and deleting of models after at a clean shut down");
				optp.println("#interactive_cloning_session\t" + interactive_cloning_session);
				if(tipps) optp.println("with blobs all neighbors of a changed pixel also change by 1/2");
				optp.println("#blobs\t" + blobs);
				if(tipps) optp.println("either resume with the last activity, or reset to fresh centered state");
				optp.println("#reset_activity\t" + reset_activity);
				if(tipps) optp.println("if input dialogs of master dataset should be invoked even during batchrun");
				optp.println("#master_shows_dialog\t" + master_shows_dialog);
				if(tipps) optp.println("perform weight optimization after every change in a model");
				optp.println("#aggressive_weights\t" + aggressive_weights);		
			}
			
			if(tipps) optp.println("create persistent backup Image of Idata before noisification");
			optp.println("#duplicate_Idata\t" + duplicate_Idata);	
			if(tipps) optp.println("scan network for available workers, disable it when preparing/launching job files");
			optp.println("#check_free_nodes\t" + check_free_nodes);	
			if(tipps) optp.println("If true the seed from this file will be used to initialize the random generator");
			optp.println("#previous_seed\t" + previous_seed);
			if(master) 
			{
				optp.println("debug_levels 0 .. none, 1 .. Pval.txt, 2 .. Autosave every wrap, 3 .. Cachreports, 4 .. Job.bin ");
				optp.println("#debug_level\t" + debug_level);
				
				optp.println("Experts/Debug section theses are not accessible from the gui");
				if(tipps) optp.println("Pixeladjustment may use quad interpolation or simple bisection");
				optp.println("#enable_quad_fits\t" + enable_quad_fits);
				if(tipps) optp.println("Pixeladjustment stops after first bisection");
				optp.println("#rough_optimization\t" + rough_optimization);
				if(tipps) optp.println("Parameters for noise in reseted or resumed models");
				optp.println("#noisify_at_restart\t" + noisify_at_restart);
				optp.println("#area_noise_level\t" + area_noise_level);
				optp.println("#noisification_level\t" + noisification_level);

				if(tipps) optp.println("the weight for a model that is Reset");
				optp.println("#default_weight\t" + default_weight);
				if(tipps) optp.println("fair weights are equal for all models except model0");
				optp.println("#fair_weights\t" + fair_weights);
			}
			if(tipps) optp.println("these control the scaling of the Ptable (VERY CRITICAL!)");
			optp.println("#pvaloffset\t" + pvaloffset);
			optp.println("#autoptableoffset\t" + autoptableoffset);
			if(master) 
			{
				optp.println("#looptotal\t" + looptotal);
				optp.println("#frame_shuffle_seed\t" + seed);
				if(tipps) optp.println("formating options for console output");
				optp.println("#report_interval\t" + report_interval);
				optp.println("#report_num\t" + report_num);
				optp.println("#report_cols\t" + report_cols);
				optp.println("#weight_cols\t" + weight_cols);
				optp.println("#max_supported_workers\t" + MAX_THREADS);
			}
			optp.flush();
		}
		catch (Exception e)
		{
			System.out.println("Error: could not write to "+ options.getPath());
			e.printStackTrace();
		}
	}

	void saveJobrun(File jobrun)
	{
		try
		{
			java.io.FileOutputStream opts = new java.io.FileOutputStream(jobrun, false);
			java.io.PrintWriter optp = new java.io.PrintWriter(opts);
			optp.println("options are read from lines that start with #");
			optp.println("These settings MUST match the prepared jobfiles for all the workers");
			optp.println("They will override any conflicting ones in options.txt");
			optp.println("job files were created from #Input_Data\t" + Datatitle);
			optp.println("located in #sourcepath\t" + sourcepath);
			optp.println("#modelnum\t" + modelnum);
			optp.println("#modelsize\t" + modelsize);
			optp.println("#modelGraylevels\t" + modelGraylevels);
			optp.println("#bondlength\t" + bondlength);
			optp.println("#data_bg0\t" + data_bg0);
			optp.println("#data_bg_std0\t" + data_bg_std0);
			optp.println("#datamean0\t" + datamean0);
			optp.println("#modelmin0\t" + modelmin0);
			optp.println("#modelmax0\t" + modelmax0);
			optp.println("#use_max\t" + use_max);
			optp.println("#use_weights\t" + use_weights);
			
			optp.flush();
		}
		catch (Exception e)
		{
			System.out.println("Error: could not write to "+ jobrun.getPath());
			e.printStackTrace();
		}
	}

	boolean loadJobrun(File jobrun)
	{
		try
		{
			if(jobrun.isFile() && jobrun.setReadable(true) )
			{
				Scanner fs = new Scanner(new FileReader(jobrun));
				while( fs.hasNextLine() )
				{
					String line = fs.nextLine();
					if(!line.equals(""))
					{
						Scanner ls = new Scanner(line);
						String token = ls.next();
						token = token.trim();
						if (token.startsWith("#"))
						{
							if(token.equals("#modelnum"))
							{	modelnum = ls.nextInt();}
							else if(token.equals("#modelsize"))
							{	modelsize = ls.nextInt();}
							else if(token.equals("#bondlength"))
							{	bondlength = ls.nextInt();}
							else if(token.equals("#modelGraylevels"))
							{	modelGraylevels = ls.nextInt();}
							else if(token.equals("#data_bg0"))
							{	data_bg0 = ls.nextDouble();}
							else if(token.equals("#data_bg_std"))
							{	data_bg0 = ls.nextDouble();}
							else if(token.equals("#datamean0"))
							{	datamean0 = ls.nextDouble();}
							else if(token.equals("#modelmin0"))
							{	modelmin0 = ls.nextDouble();}
							else if(token.equals("#modelmax0"))
							{	modelmax0 = ls.nextDouble();}
							else if(token.equals("#use_weights"))
							{	
								boolean b = ls.nextBoolean();
								tell_use_weights = (b != use_weights);	
							}	
							else if(token.equals("#use_max"))
							{	
								boolean b = ls.nextBoolean();
								tell_use_max = (b != use_max);
							}		
						}
						ls.close();
					}
				}
				fs.close();
			}
		}
		catch (Exception e)
		{
			System.out.println("Error: could not read " + jobrun.getPath());
			e.printStackTrace();
			return false;
		}
		return true;
	}

	void saveWeights(String my_path)
	{
		/*
		int mnum = modelnum;
		if(cloned)
		{	mnum = (modelnum+1)/2;} //either doubled or doubled - 1 
		*/
		try
		{
			java.io.FileOutputStream weightStream = new java.io.FileOutputStream(my_path + "Weights.txt", false);
			java.io.PrintWriter weightWriter = new java.io.PrintWriter(weightStream);
			weightWriter.println("Any line starting with an Integer will be interpreted");
			weightWriter.println("as a modelnumber, a modelweight and a Modelstatus");
			weightWriter.println("Beware, the modelnumbers themselfes are actually ignored, only the valid line count matters");
			weightWriter.println("Lattice\t\timmutable model with fixed weight");
			weightWriter.println("Continue or \"blank\"\tmutable model with adjustable weight");
			weightWriter.println("Finished\t\tconverged model with adjustable weight");
			weightWriter.println("Fixed\t\tmutable model with fixed weight");
			weightWriter.println("Reset\t\tA new model as a noisy clone of the Lattice and then Continued");
			weightWriter.println("Noisify\t\tNoise is added to the model and then it is Continued");
			weightWriter.println("After modifications you have to save and tick reload weights");
			for (int m=0; m < modelnum; ++m)
			{
				/*
				double weight = modelweight[i];
				double scaling = cloned ? (1.0-clone_weight) : 1.0;
				if( (i==0) && skip_first_clone && (mnum>1) ) //we could but that into the line above
				{	scaling = 1.0;} 
				weightWriter.println("" + (i+1) + "\t" + ( scaling*weight ) + "\t" + (cloned?"Continue":getModelStateString(dead[i])));
				*/
				weightWriter.println("" + (m+1) + "\t" + modelweight[m] + "\t" + getModelStateString(dead[m]));
				//System.out.println("" + (m+1) + "\t" + modelweight[m] + "\t" + getModelStateString(dead[m]));
			}
			/*
			if( cloned )
			{
				
				int m0 = (( skip_first_clone&&(mnum>1) )?1:0);
				for (int i = m0; i < mnum; ++i)
				{
					double weight = modelweight[i];
					double scaling = clone_weight;
					weightWriter.println("" + (i+1+mnum-m0) + "\t" + (scaling*weight) + "\t" + "Continue");
					//System.out.println("" + (i+1) + "\t" + modelweight[i] + "\t" + getModelStateString(dead[i]));
				}
			}
			*/
			weightWriter.close();
			weightStream.close();
			Weighttitle = "Weights.txt"; //from now on we would rather read the file
			System.out.println("saved Weights.txt"); 
		}
		catch (Exception e)
		{
			System.out.println("Error: could not write to " + my_path + "Weights.txt");
			e.printStackTrace();
		}
	}

	void saveLog(String path)
	{
		try
		{
			java.io.FileOutputStream logStream = new java.io.FileOutputStream(path + "Log.txt", true);
			java.io.PrintWriter logWriter = new java.io.PrintWriter(logStream);
			logWriter.println("******* NEW SESSION *******");
			String log = IJ.getLog();
			logWriter.write(log);
			logWriter.close();
			logStream.close();
		}
		catch (Exception e)
		{
			System.out.println("Error: could not write to " + path + "Log.txt");
			e.printStackTrace();
		}
	}

	boolean writeJobStats(File jobStats)
	{
		try
		{
			DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
			Date date = new Date();
			date.setTime(System.currentTimeMillis());
			
			java.io.FileOutputStream jobStatStream = new java.io.FileOutputStream(jobStats, false);
			java.io.PrintWriter jobStatWriter = new java.io.PrintWriter(jobStatStream);
			
			jobStatWriter.println( "created by " + user + " at " + dateFormat.format(date) );
			jobStatWriter.println("#numjobs\t" + 1); //jobVec.size()
			//for(int job = 0; job < jobVec.size(); ++job)
			//{
				double modelmin = jobVec.get(default_job).modelmin;
				double modelmax = jobVec.get(default_job).modelmax;
				long datasum = jobVec.get(default_job).datasum;
				int dataStSize = jobVec.get(default_job).use_frames;
				//long[] dataSumImg = jobVec.get(job).dataSumImg;
				
				jobStatWriter.println("#job\t" + default_job);
				jobStatWriter.println("#modelmin\t" + modelmin);
				jobStatWriter.println("#modelmax\t" + modelmax);
				jobStatWriter.println("#datasum\t" + datasum);
				jobStatWriter.println("#dataStSize\t" + dataStSize);
				jobStatWriter.println("#total_subframes\t" + total_subframes);
				jobStatWriter.println("#total_SumImg\t" + bigSumImg.length);
				for(int pix = 0; pix < bigSumImg.length; ++ pix)
				{
					jobStatWriter.print("" + bigSumImg[pix]);
					if((pix+1) % 20 != 0)
					{
						jobStatWriter.print("\t");
					}
					else
					{	
						jobStatWriter.println();
					}
				}
				jobStatWriter.println();
				
			//}
			jobStatWriter.close();
			jobStatStream.close(); 	
		}
		catch (Exception e)
		{
			System.out.println("Error: could not write to " + jobStats.getPath());
			e.printStackTrace();
			return false;
		}
		return true;
	}

	boolean readJobStats(File jobStats)
	{
		if(! (jobStats.isFile() && jobStats.setReadable(true) ) )
		{	return false;}
		try
		{
			Scanner fs = new Scanner(new FileReader(jobStats));
			while( fs.hasNextLine() )
			{
				//FIXME Why does reading job from file fail?
				//int job = 0; //to be read from file //crashes with -1!?
				int num_pix = -1; //has to be read from file
				String line = fs.nextLine();
				if(!line.equals(""))
				{
					Scanner ls = new Scanner(line);
					String token = ls.hasNext()?ls.next():"error";
					//String token = ls.next();
					token = token.trim();
					if (token.startsWith("#"))
					{
						if(token.equals("#numjobs"))
						{	
							int numjobs = ls.nextInt();
							while(jobVec.size() < numjobs)
							{	
								jobVec.add( new JobData() ); 
							}
						}
						else if(token.equals("#job"))
						{	
							default_job = ls.nextInt();
						}
						else if(token.equals("#modelmin"))
						{	
							jobVec.get(default_job).modelmin = ls.nextDouble();
						}
						else if(token.equals("#modelmax"))
						{	
							jobVec.get(default_job).modelmax = ls.nextDouble();
						}
						else if(token.equals("#datasum"))
						{	
							datasum0 = ls.nextLong();
							jobVec.get(default_job).datasum = datasum0;
						}
						else if(token.equals("#total_subframes"))
						{	
							total_subframes = ls.nextInt();
							IJ.log("total frames: " + total_subframes);
						}	
						else if(token.equals("#dataStSize"))
						{	
							dataStSize0 = ls.nextInt();
							jobVec.get(default_job).use_frames = dataStSize0;
						}	
						else if(token.equals("#dataSumImg") || token.equals("#total_SumImg"))
						{	
							num_pix = ls.nextInt();
							jobVec.get(default_job).dataSumImg = new long[num_pix];
							bigSumImg = new long [num_pix];
							for(int ind=0; ind < num_pix; ++ind)
							{
								long val = fs.nextInt();
								jobVec.get(default_job).dataSumImg[ind] = val; //There may be several per line
								bigSumImg[ind] = val;
							}	
						}					
					}
					ls.close();
					if(total_subframes < 1)
					{
						total_subframes = dataStSize0; //correct fix for old-style single job data
					}
				}
			}
			fs.close();
		}
		catch (Exception e)
		{
			System.out.println("Error: could not read " + jobStats.getPath());
			e.printStackTrace();
			return false;
		}
		return true;
	}


	void writeSignificance(String path, ImagePlus Models)
	{
		try
		{
			DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
			Date date = new Date();
			date.setTime(System.currentTimeMillis());
			java.io.FileOutputStream sigStream = new java.io.FileOutputStream(path + "Significance.txt", true);
			java.io.PrintWriter sigWriter = new java.io.PrintWriter(sigStream);
			sigWriter.println( dateFormat.format(date) );
			ImageStack ModelsSt = Models.getStack();
			for(int job = 0; job < jobVec.size(); ++job)
			{
				double modelmin = jobVec.get(job).modelmin;
				double modelmax = jobVec.get(job).modelmax;
				long datasum = jobVec.get(job).datasum;
				int dataStSize = jobVec.get(job).use_frames;
				long[] dataSumImg = jobVec.get(job).dataSumImg;
				
				sigWriter.println("");
				sigWriter.println("job" + job + "\t" + ( (jobVec.get(job).Idata==null)?"":jobVec.get(job).Idata.getTitle() ) + "\tframes: " + dataStSize);
				sigWriter.println( "total Data Counts: " + datasum);
				
				int modelcounts = 0;
				int[] devAvg = new int[modelnum];
				int[] devPlus = new int[modelnum];
				int[] devMinus = new int[modelnum];
				int[] counts = new int[modelnum];
				double[] cts_case = new double[modelnum];
				double[] delta_p_case = new double[modelnum];
				double[] delta_m_case = new double[modelnum];
				for(int m = 0; m < modelnum; ++m)
				{
					short[] pixels = (short[])ModelsSt.getPixels(1+m);
					int pixsum = 0;
					double devSum = 0.0;
					double plusSum = 0.0;
					double minusSum = 0.0;
					double wght = modelweight[m];
					for(int i = 0; i < pixels.length; ++i )
					{
						pixsum += (int)pixels[i];
						double signal = (modelmin + (double)pixels[i] * modelmax / modelGraylevels) * dataStSize ;
						devSum += Math.abs( signal - dataSumImg[i] );
						if(signal-dataSumImg[i] > 0)
						{
							plusSum += (signal-dataSumImg[i]);
						}
						else
						{
							minusSum += (dataSumImg[i] - signal);
						}
					}
					cts_case[m]     = devSum   / (double)dataStSize; 
					delta_p_case[m] = plusSum  / (double)dataStSize;
					delta_m_case[m] = minusSum / (double)dataStSize; 
					
					devAvg[m]   = (int) Math.round(devSum   * wght);
					devPlus[m]  = (int) Math.round(plusSum  * wght);
					devMinus[m] = (int) Math.round(minusSum * wght);
					counts[m] = (int) Math.round( (modelmin + modelmax * ((double)pixsum) / modelGraylevels) * wght * dataStSize);
					modelcounts += counts[m];		
				}
				
				sigWriter.println("total Model Counts: " + modelcounts);
				sigWriter.println("");
				for(int m = 0; m < modelnum; ++m)
				{
					double wght = modelweight[m];
					double rel_signal = (double) devAvg[m] / (double) modelcounts;
					
					int cases2 = (int) Math.round( dataStSize * wght );
					sigWriter.println("Model: " + m + "\tWeight: " + format.format(wght) +
									  "\tCounts: " + counts[m] + "\t\tD Counts: " + devAvg[m] );
					sigWriter.println("Model: " + m + "\tDabs: " + format.format(rel_signal) +
					"\t\tcases: " + cases2 + "\t\tDabs/case: " + format.format(cts_case[m]) );
					sigWriter.println("Model: " + m + "\tD+: " + devPlus[m] +
					"\t\tD-: " + devMinus[m]);
					sigWriter.println("Model: " + m + "\tD+/case: " + format.format(delta_p_case[m]) +
					"\tD-/case: " + format.format(delta_m_case[m]) + "\tD/case: " + format.format(delta_p_case[m]-delta_m_case[m]));
					
					
					if(m != modelnum-1)
					{	sigWriter.println("");}
				}
				sigWriter.println("");
			}
			sigWriter.close();
			sigStream.close();
		}
		catch (Exception e)
		{
			System.out.println("Error: could not write to " + path + "Significance.txt");
			e.printStackTrace();
		}
	}

	void writePerformance(String path)
	{
		try
		{
			java.io.FileOutputStream perfStream = new java.io.FileOutputStream(path + "Performance.txt", false);
			java.io.PrintWriter perfWriter = new java.io.PrintWriter(perfStream);
			perfWriter.println("multiplicity\tperformance\tlocal\texecutable");
			for(int i = 0; i < use_threads; ++i)
			{
				int job = threadVec.get(i).job;
				int rank = threadVec.get(i).rank;
				double perf = jobVec.get(job).performances[rank];
				boolean local = threadVec.get(i).can_write;
				
				String exec = threadVec.get(i).command;
				perfWriter.println("1\t" + perf + "\t" + local + "\t" + exec);
				perfWriter.flush();
				perfStream.flush();
			}
			perfWriter.close();
			perfStream.close();
		}
		catch (Exception e)
		{
			System.out.println("Error: could not write to " + path + "Performance.txt");
			e.printStackTrace();
		}
	}
	
	void sync_images(ImagePlus sync_imps[])
	{
		ImagePlus imp = WindowManager.getCurrentImage();		
		if (imp != null)
		{
			int slice = imp.getSlice();
			if(slice != the_slice)
			{
				boolean on_list = false;
				for(int i = 0; i < sync_imps.length; ++i)
				{
					if(imp==sync_imps[i])
					{
						on_list = true;
						break;
					}
				}
				if(on_list)
				{
					the_slice = slice;
					for(int i=0; i < sync_imps.length; ++i)
					{
						sync_imps[i].setSlice(the_slice);
					}
				}
			}
		}
	}

	boolean parse_batchfile(int lvl, File batchfile)
	{	
		final int max_lvl = 3; //limit for nested includes
		//IJ.log("parsing (level: " + lvl + "): " + batchfile.getPath() );
		if(lvl == 0)
		{
			IJ.open(batchfile.getPath());
			IJ.selectWindow(batchfile.getName());
			NonBlockingGenericDialog bdg = new NonBlockingGenericDialog(getClass().getSimpleName() + " batchmode");
			bdg.addMessage("You may now edit " + batchfile.getName() + " and save your changes");
			bdg.addMessage("Click ok to re-read " + batchfile.getName() + " from disk and proceed");
			bdg.addMessage("The batchrun wont continue if " + batchfile.getName() + " is closed");
			bdg.addMessage("You may halt and modify a running reconstruction with Esc");
			
			bdg.showDialog();
			
			if(bdg.wasCanceled())
			{	return false;}
		}
		try
		{
			Scanner fs = new Scanner(new FileReader(batchfile));
			int repeat = 1;
			boolean skipping = false;
			boolean syntax_error = false;
			while( fs.hasNextLine() && !skipping)
			{
				String line = fs.nextLine();
				line = line.trim();
				if(!line.equals("") && !line.startsWith("#"))
				{
					boolean is_a_path = true;
					//boolean is_a_file = false;
					if(line.startsWith("&"))
					{
						is_a_path = false;
						Scanner ls = new Scanner(line);
						String command = ls.next().trim();
						//IJ.log("command: " + command);
						if(command.equals("&REPEAT"))
						{
							if(ls.hasNextInt())
							{	repeat = ls.nextInt();}
							else
							{	
								IJ.log("Syntax Error in " + batchfile.getPath() + " Integer expected after &REPEAT");
								syntax_error = true;	
							}
							//is_a_path = false;	
						}
						else if(command.equals("&PREFIX"))
						{
							if(ls.hasNextLine())
							{	prefix = ls.nextLine().trim();}
							else
							{	prefix = "";}
							//is_a_path = false;
						}
						else if(command.equals("&OPTIONS")) //syntactically equivalent to naked filename
						{
							boolean line_ok = false;
							if(ls.hasNextLine())
							{	
								line = ls.nextLine().trim();
								if(!line.equals(""))
								{	line_ok = true;}
								if(!line_ok)
								{	
									IJ.log("Syntax Error in " + batchfile.getPath() + " string (filename) expected after &OPTIONS");
									is_a_path = false;
									syntax_error = true;	
								}
									
							}
							
							is_a_path = true;
							//is_a_file = true;
						}
						else if(command.equals("&HERE"))
						{	
							line = batchfile.getParent() + File.separator;
							if(ls.hasNextLine())
							{	line += ls.nextLine().trim();}
							//IJ.log("HERE: " + line);
							is_a_path = true;
						}
						else if(command.equals("&EXIT"))
						{
							skipping = true;
							//is_a_path = false;
						}
						else if(command.equals("&INCLUDE"))
						{
							//is_a_path = false;
							boolean line_ok = false;
							
							line = ( ls.hasNextLine() ? ls.nextLine().trim() : "");
							if(!line.equals(""))
							{	line_ok = true;}
							if(!line_ok)
							{	
								IJ.log("Syntax Error in " + batchfile.getPath() + " string (filename) expected after &INCLUDE");
								syntax_error = true;	
							}
							else
							{
								File batchfile2 = new File(line);
								if(batchfile2.isFile() && batchfile2.setReadable(true) && (lvl < max_lvl))
								{
									for(int i=0;i<repeat;++i)
									{
										boolean include_ok = parse_batchfile(lvl + 1, batchfile2);
										if(!include_ok)
										{	return false;}
									}
									repeat = 1;
								}
								else
								{	IJ.log("Failed/refused to include " + batchfile2.getPath() + " at include level: " + (lvl+1) + " of " + max_lvl);}
							}
						}
						else if( command.equals("&OUTPUT"))
						{
							//is_a_path = false;
							if( lvl == 0)
							{
								if(ls.hasNextLine())
								{
									outpath = ls.nextLine().trim();
									if(outpath.equals("&HERE"))
									{	
										outpath = batchfile.getParent() + File.separator;
									}
									File checkpath = new File(outpath.substring(0,outpath.length()-1));
									//IJ.log("checking: " + checkpath.getPath());
									if(!checkpath.isDirectory() || !checkpath.canWrite())
									{
										IJ.log("invalid outpath: " + checkpath.getPath() + File.separator);
										outpath = "";
										syntax_error = true;
									}
								}
							}
							else
							{
								IJ.log("ignoring command &OUTPATH at include level " + lvl);
							}
						}
						else
						{	
							IJ.log("Syntax Error in " + batchfile.getPath() + " unknown command: " + command);
							//is_a_path = false;
						}
					}
					if(is_a_path)
					{
						if(line.endsWith(File.separator) )
						{	line += "options.txt";}
						for(int i=0;i<repeat;++i)
						{	
							batchpaths.add(line);
							prefix_per_path.add(prefix);
							//IJ.log("adding job (level: " + lvl + "): " + line + " prefix: " + prefix);
						}
						repeat = 1;
					}
				}
			}
			if(syntax_error)
			{	return false;}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			IJ.error("Error reading " + batchfile.getPath());
			return false;
		}		
		boolean paths_ok = true;
		if(lvl == 0) //only the top level beatchrun check all files once in the end
		{
			for(int batch_id = 0; batch_id < batchpaths.size(); ++batch_id)
			{
				File peek = new File(batchpaths.get(batch_id));
				if(!peek.isFile() && peek.setReadable(true) )
				{
					IJ.log("missing file: " + peek.getPath());
					paths_ok = false;
				}
				else
				{
					IJ.log("confirmed: " + peek.getPath());
				}
			}
		}
		return paths_ok;
	}
	
	private int switch_pixels( ImagePlus Models, int[] chx, int[] chy, int[] chm, int[] newv, int[] oldv)
	{
		int maxlen = chx.length;
		int bx = (homex + blockx + modelsize) % modelsize;
		int by = (homey + blocky + modelsize) % modelsize;
		int mm = (homem+blockm) % modelnum;
		
		ImageStack ModelsSt = Models.getStack();
		short[] pixels = (short[])ModelsSt.getPixels( mm + 1 );
	
		
		int x = bx - (by-(by&1))/2;
		int z = by;
		int y = -x-z;
		int[] cube = {x-cx,y-cy,z-cz};
		hexagonaltrap(cube, modelsize/2);
		x = cube[0] + cx;
		y = cube[1] + cy;
		z = cube[2] + cz;
		int oq = x + (z-(z&1))/2;
		int or = z;
		
		int oind = oq+or*modelsize;
		int nq = oq;
		int nr = or;
		int nind = oind;
		
		int bl = bondlength;
		int watchdog = 10;
		
		do
		{
		
			int dx = (random.nextBoolean()?1:-1) * random.nextInt(bl+1);
			int dy = (random.nextBoolean()?1:-1) * (random.nextInt(dx==0?bl:bl+1)+(dx==0?1:0));
			int dz = -dx-dy;
			int nx = x + dx;
			int nz = z + dz;
			int ny = -nx -nz;
			int[] ncube = {nx-cx,ny-cy,nz-cz};
			hexagonaltrap(ncube, modelsize/2);
			nx = ncube[0] + cx;
			ny = ncube[1] + cy;
			nz = ncube[2] + cz;
			nq = nx + (nz-(nz&1))/2;
			nr = nz;
			nind = nq + nr * modelsize;
			if (--watchdog < 0)
			{	return 0;}
		}
		while(pixels[oind]==pixels[nind]);
		//suggest to swap intensities of pixel 1 and 2
		
		chx[0] = oq;
		chy[0] = or;
		chm[0] = mm;
		oldv[0] = pixels[oind];
		newv[0] = pixels[nind];
		if(maxlen < 2) //just in case
		{	return 1;}
		
		chx[1] = nq;
		chy[1] = nr;
		chm[1] = mm;
		oldv[1] = pixels[nind];
		newv[1] = pixels[oind];
		
		return 2;
	}
	
	private int blobify_changes( ImagePlus Models, int chlistlen, int[] chx, int[] chy, int[] chm, int[] newv, int[] oldv)
	{
		
		int maxlen = chx.length;
		int[] chx2 = new int[maxlen];
		int[] chy2 = new int[maxlen];
		int[] chm2 = new int[maxlen];
		int[] newv2 = new int[maxlen];
		int[] oldv2 = new int[maxlen];
		
		
		int q2 = 0;
		//always add original points to front of blobbed changes
		for (int q = 0; q < chlistlen; ++q)
		{
			int mm = chm[q];
			short[] pixels= (short[])Models.getStack().getPixels(mm+1);
			short ov = pixels[ chx[q] + modelsize*chy[q] ];
			short nv = (short)newv[q];
			if(q2 < maxlen)
			{
				chx2[q2] = chx[q];
				chy2[q2] = chy[q];
				chm2[q2] = chm[q];
				newv2[q2] = newv[q];
				oldv2[q2] = ov;
			}
			++q2;
		}  
		
		for (int q = 0; q < chlistlen; ++q)
		{
			int mm = chm[q];
			short[] pixels= (short[])Models.getStack().getPixels(mm+1);
			short ov = pixels[ chx[q] + modelsize*chy[q] ];
			short nv = (short)newv[q];
			//int bind = blockx + blocky * blocksize + blockm * blockarea;
			if(blobs && ( (nv-ov)/2 != 0 ) )
			{
				int mq = chx[q];
				int mr = chy[q];
				
				int x = mq - mr/2;
				int z = mr;
				int y = -x -z;
				
				for(int dx = -1; dx <= 1; ++dx)
				{
					for(int dz = -1; dz <= 1; ++dz)
					{
						int dy = -dx - dz;
						if(Math.abs(dx)+Math.abs(dy)+Math.abs(dz) == 2)
						{
							int nx = x + dx;
							int nz = z + dz;
							int ny = -nx -nz;
							int[] cube = {nx-cx,ny-cy,nz-cz};
							hexagonaltrap(cube, modelsize/2);
							nx = cube[0] + cx;
							ny = cube[1] + cy;
							nz = cube[2] + cz;
							int nq = nx + nz/2;
							int nr = nz;
							int ind = nq + nr * modelsize;
							short oNv = pixels[ ind ];
							int bind = blockx + blocky * blocksize + blockm * blockarea;
							/**   only round at initial guess and simple step search**/
							short nNv = (short)( (incrs[bind]<2)  ? 
												 ((oNv+nv)/2) :
												 (oNv + (nv-ov)/2)	);
							//short nNv = (short)((oNv+nv)/2); /**  always rounding**/
							
							//short nNv =(short) (((oNv+nv)+1)/2);
							if(nNv <= 0)
								{nNv = 0;}
							if(nNv >= ( modelGraylevels-1 ) )
								{nNv = (short) (modelGraylevels-1);}
							//add the Neigboring point to the blobbed list
							if(nNv != oNv)
							{
								int qf = -1; //we dont know that place yet
								for(int qs=0;(qs<q2) && (qs < maxlen);++qs) 
								{
									if( (nq==chx2[qs]) && (nr==chy2[qs]) && (mm==chm2[qs]) )
									{
										qf = qs;
										break;
									}	  
								} 
								
								
								
								if(qf == -1 ) // add the point
								{								
									if(q2 < maxlen)
									{
										chx2[q2] = nq;
										chy2[q2] = nr;
										chm2[q2] = mm;
										newv2[q2] = nNv;
										oldv2[q2] = oNv;
									}
									++q2;
								}
								else //update it
								{
									newv2[qf] = newv2[qf] + nNv - oNv;
								}
							}  
						}	
					}
				}
			}
		}
		///WARNING there is no check for duplicated neighbors //should work now
		
		
		int q3 = 0;
		for(int q=0; q<q2;++q)
		{
			if(newv2[q]!=oldv2[q])
			{
				if(q3 < maxlen)
				{
					chx[q3] = chx2[q];
					chy[q3] = chy2[q];
					chm[q3] = chm2[q];
					newv[q3] = newv2[q];
					oldv[q3] = oldv2[q];
				}
				++q3;
			}
		}
		return q3; //the new chlistlen
	}
	
	
	private double make_nextModel(ImagePlus Models, int chlistlen, int[] chx, int[] chy, int[] chm, int[] newv, int[] oldv)
	{
		///WARNING no check that all changes are actually on the same model
		int mm = chm[0];
		short[] pixels= (short[])Models.getStack().getPixels(mm+1);
		short[] nextpixels= (short[])nextModel.getStack().getPixels(1);
		for(int i = 0; i< pixels.length; ++i)
		{
			nextpixels[i] = pixels[i];
		}
		
		for(int q = 0; q < chlistlen; ++q)
		{
			int ind = chx[q] + modelsize * chy[q];
			nextpixels[ind] = (short)newv[q];
		}
		if(nextMatcher == null)
		{
			return 1.0;
		}
		else
		{
			nextMatcher.update_matches();
			return median_filtering ? nextMatcher.summed_matches[0] : nextMatcher.inv_avg_matches[0];
		}
	}
	
	boolean do_Model_Manip_Dialog( int[] vals)
	{
		boolean first = true;
		do
		{
			NonBlockingGenericDialog gd = new NonBlockingGenericDialog(getClass().getSimpleName() + " Model Manipulation");
			gd.setSmartRecording(true);
			gd.addMessage("Click Ok to commit the operations\n"+
			"Click Cancel once you are done");
			gd.addNumericField("first model(1.." + modelnum + ")", vals[0], 0);
			gd.addNumericField("last model(1.." + modelnum + ")", vals[1], 0);
			gd.addNumericField("children( 0..delete, 1..swap, 2..clone )", vals[2], 0);
			gd.addNumericField("clone weight", clone_weight, 4);
			if(!first)
			gd.addMessage("Please Check Your Inputs");
			gd.showDialog();
			if(gd.wasCanceled())
			{
				return false;
			}
			vals[0] = (int)gd.getNextNumber();
			vals[1] = (int)gd.getNextNumber();
			vals[2] = (int)gd.getNextNumber();
			clone_weight = gd.getNextNumber();
			first = false;
		}
		while(!( (clone_weight>=0.0) && 
				 (vals[0] > 0) && (vals[0] <= modelnum) &&
				 (vals[1] > 0) && (vals[1] <= modelnum) &&
				 (vals[2] >=0) ));
		return true;
	}
	
	void modelManipulation(boolean force_dialog)
	{
		int fi = skip_first_clone?2:1;
		int la = modelnum;
		int num_cl = 2; //0..delete, 1..no operation, 2 .. 2 copies, ...
		boolean imm = interactive_cloning_session || force_dialog;
		do
		{
			int pars[] = {fi,la,num_cl};
			if(imm)
			{
				if( !do_Model_Manip_Dialog(pars) )
				{	break; }
			}
			fi = pars[0];
			la = pars[1];
			num_cl = pars[2];
			int[] nxt_seq = nextModelSequence(fi, la, num_cl);
			
			if( nxt_seq != null )
			{
				Activity = expandImg( Activity, nxt_seq );
				ActivitySt = Activity.getStack();
				Updates = expandImg( Updates, nxt_seq );
				UpdatesSt = Updates.getStack();
				Models = expandImg( Models, nxt_seq );
				ModelsSt = Models.getStack();
				Histimg = expandImg( Histimg, nxt_seq );
				//this stack is rarely used and always local
				Diffimg = expandImg( Diffimg, nxt_seq );
				DiffSt = Diffimg.getStack();
				expandWeights( fi, la, num_cl);
			}	
			IJ.log("new modelnum: " + modelnum);
		}
		while(imm);
	}
	
}
