import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.io.FileInfo;
import ij.text.TextWindow;
import ij.plugin.*;
import java.awt.*;
import java.awt.geom.*;
import java.awt.Rectangle;

import java.util.Random;
import java.util.Scanner;
import java.util.Vector;
import java.util.Locale;
import java.util.HashMap;
import java.util.Properties;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.PrintWriter;
import java.io.IOException;

import java.nio.ByteBuffer;
import java.text.DecimalFormat;

import java.util.Date;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

public class Hex_Sampler implements PlugIn 
{
	static Random random = null;
	static byte[] bytes = new byte[8]; 
	static ByteBuffer buffer = ByteBuffer.wrap(bytes);
	static java.text.DecimalFormat format = new DecimalFormat("0.000000");
    static DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	static Date date = new Date();
	static int mutex = 0;
	
	static byte flag_hole = 0;
	static byte flag_graphene = 1;		
	static byte flag_molecule = 4;
	static byte flag_dirt = 16;
	static byte flag_taken = 64;

    
    String[] filenames = {"frame_init.txt","frame_continue.txt","frame_sorted.csv","<open>","<use defaults>"}; 
	String[] filter_choices = {"<none>","<open>","filter_frames.txt"};
	String[] op_modes = {"optimize merit","cut subframes","poll subframes","debug run","quick"};
    //String[] merits = {"default c * m * l","hexagonal (c in hexpixels)"};
    String[] states = {"skipped","done","pending","init", "busy", "done/redo"};
    String prefix = "";
    Vector<StringBuilder> groupHeaders = new Vector<StringBuilder>();
    
    int mask_scaling = 1;
    int stability = 1;
    int peak_points = 5;
    int initthreadnum = 0;
    int total_subframes = 0;
    int threadnum = initthreadnum;
    int streamnum = 0;
    int thread_offset = 0; //typically hexclient would busy the first 2 slots per PC
    int bondlength = 4; //hexpixels per bondlength
    int minilength = 12;
    int max_threads = 0; //soft limit
    int target_updates = 0;
    int deaths = 0;
    int max_relaunches = 0; //abandon a worker after 2 exceptions/crashes/connection losses whatever error
	
	boolean retune_origin = false;
	boolean locking = false; //if workes on the same computer should enque for responding 
	double min_molecule_coverage = 0.0;
	double max_molecule_coverage = 0.0;
	int molecule_value = 4;
	
	int max_stack_size = 1000000;
	int stack_counter = 0;
	int init_use_frames = 0;
	int use_frames = init_use_frames;
	int frames_done = 0;
	int uc_res = 32;
	int smooth_passes = 1;
	int monitor_worker = 0;
	int op_mode = 0; // 0 .. optimize merit, 1 .. cut subframes, 2 .. poll subframes 3 .. debug 4 .. quick
	int the_slice = 1;
	int total_resets = 0;
	
	int next_frame = 0;
	//int finished_frame = 0;
	int active_frames = 0;
	//int[] busy_reports = null;
	int[] busy_frames = null; //should_become cur_rep
	int[] fr_p_worker = null;
	int modelsize = 64;
	int solidsize = modelsize/2;
	double field_of_view = 20;
	final double a0 = 0.142;
	double hexedgelen = 1.66; //typical value
	double def_hel = 0.0;
	double def_excent = 1.03;
	double def_phi = 0.0;
	double min_hex_contrast = 0.01;
	
	double tilt = 0.0;
	double mean_tilt = 0.0;
	double mean_hel = 0.0;
	double mean_ellA = 0.0;
	double mean_ellB = 0.0;
	double mean_phi = 0.0;
	boolean keep_Idata = true;
	boolean check_free_nodes = false;
	boolean forever = false; //cant be loaded from options.txt as it defies the purpose of a batchrun
	boolean endless = false;
	boolean reorder_stack = false;
	boolean backwards = false;
	boolean use_mirrorQ = true;
	boolean use_positionQ = true;
	boolean binary_mode = true; //much quicker at launching but also binary job files 
	boolean reset_states = true;
	boolean mark_sf = false;
	boolean reset_bad_frames = true;
	boolean has_mask = false;
	boolean skip_grid = true;
	int create_defect = 0;
	int rotate_tilt = 0;
	boolean new_submissions = true;
	int debug_level = 0;
	double total_delta_merit = 0.0;
	double total_delta_hex_merit = 0.0;
	String path = IJ.getDirectory("imagej");
	String sourcepath = path; //location of the imagestacks
	
	String user = System.getProperty("user.name");
	Vector<String> batchpaths = new Vector<String>();
	File batchfile = null;
	File nodes = null;
	File options = null;
	String inputname = "frame_continue.txt";
	File frame_file = null;
	String filter_name = "<none>";
	File filter_file = null;
	//java.io.FileOutputStream resultstream = null;
	//java.io.PrintWriter resultwriter = null;
	TextWindow statusdisplay = null;
	boolean batch_run = false;
	boolean worker_ready = false;
	boolean target_useable = true;
	boolean target_needs_update = false;
	boolean update_target = false;
	
	String oldDatatitle = null;
	String Datatitle = null;
	ImagePlus Idata = null;
	ImageStack IdataSt = null;
	
	String oldMasktitle = null;
	String Masktitle = null;
	ImagePlus Mask = null;
	ImageStack MaskSt = null;
	
	ImagePlus uc_sum = null;
	ImageStack uc_sumSt = null;
	
	ImagePlus uc_sym = null;
	ImageStack uc_symSt = null;
	
	ImagePlus hexImg = null;
	ImageStack hexImgSt = null;
	
	ImagePlus sf_sum = null;
	ImageStack sf_sumSt = null;
	
	ImagePlus[] sfImgArray = null;
	ImagePlus sfImg = null;
	ImageStack sfImgSt = null;
	
	ImagePlus statImg = null;
	ImageStack statImgSt = null;
	
	ImagePlus uc_mini = null;
	ImageStack uc_miniSt = null;
	
	String uc_target_title = "<new>";
	String old_uc_target_title = "<new>";
	ImagePlus uc_target = null;
	ImageStack uc_targetSt = null;
	
	
	FrameFilter filter = null;
	
	int impHeight = 1;
	int impWidth = 1;
	int impDepth = 1;
	int mskDepth = 1;
	double hexImgScale = 1.0;
	int threadadd = 0;
	int use_threads = 0;
	Vector<ThreadData> threadVec = new Vector<ThreadData>();

	java.lang.Runtime rt = null;
	java.lang.Process[] p = null;
	java.io.BufferedInputStream[] is = null;
	java.io.DataInputStream[] dis = null;
	java.io.InputStream[] err_s = null;
	java.io.InputStreamReader[] reader = null;
	java.io.BufferedReader[] err_reader = null;
	java.io.OutputStream[] os = null;
	java.io.PrintWriter[] pwos = null;

	HexSampleReport[] total_reports = null;
	
	private void resetGlobals()
	{
		oldDatatitle = null;
		oldMasktitle = null;
		Idata = null;
		IdataSt = null;
		Mask = null;
		MaskSt = null;
		uc_target = null;
		uc_targetSt = null;
		p = null;
		is = null;
		dis = null;
		err_s = null;
		reader = null;
		err_reader = null;
		os = null;
		pwos = null;
		busy_frames = null;
		
		fr_p_worker = null;
		frame_file = null;
		filter_file = null;
		inputname = "open";
		threadVec.clear();
		total_subframes = 0;
		threadadd = 0;
		use_threads = 0;
		streamnum = 0;
		next_frame = 0;
		frames_done = 0;
		reset_states = true;
		molecule_value = 4;
		retune_origin = false;
		create_defect = 0;
		min_molecule_coverage = 0.0;
		max_molecule_coverage = 0.0;
		total_resets = 0;
		//finished_frame = 0;
		
	}
	
	
	public void run(String arg)
    {
		if(mutex == 1)
		{	
			IJ.error("Sorry, there is aleady another instance of " + getClass().getSimpleName() + " running!");
			return;
		}
		Locale.setDefault(Locale.UK);
		rt = java.lang.Runtime.getRuntime();
		
		date.setTime(System.currentTimeMillis());
        IJ.log(getClass().getSimpleName() + " " + dateFormat.format(date) + "  " + user);
		path = IJ.getDirectory(getClass().getSimpleName());
		if(path == null) {	return;}
        IJ.log(path);
        
        
        batchfile = new File(path + "batchrun.txt");
        if(batchfile.isFile() && batchfile.setReadable(true))
        {	if(!readBatchFile()) return;}   
        if(batchpaths.isEmpty())
		{	batchpaths.add(path);}
		options = new File(batchpaths.get(0) + "options.txt");
        if( !loadOptions())
        {	return;}
		
		for(int batch_id = 0; batch_id < batchpaths.size(); ++batch_id)
		{
			if(batch_id != 0)
			{
				resetGlobals();
				options = new File(batchpaths.get(batch_id) + "options.txt");
				if( !loadOptions())
				{	return;}
			}
			if(batch_run)
			{
				statusdisplay.append(options.getPath() + "(" + ( batch_id + 1 ) + " of " + batchpaths.size() + ")" );
				System.out.println("Running job " + ( batch_id + 1 ) + " of " + batchpaths.size() );
				IJ.log("Running job " + (batch_id + 1) + " of " + batchpaths.size());
			}
			IJ.log(options.getPath());
			max_threads = loadNodes();
			IJ.log("available workers: " + max_threads);
			
			if( !(  fetch_Idata() ) ) {return;}
			if( !(  fetch_Mask() ) ) 
			{
				has_mask = false;
				Masktitle = "<none>";	
			}
			if( !fetch_uc_target() )
			{	uc_target_title = "<new>";}
			int res = 0;
			do
			{
				do
				{
					if(!setupDialog())	{	return;}
					res = validateInput();
					if(res > 0)	{	IJ.error("bad input! Better luck next time");}
				} while (res > 0);
			}
			while( !initialize_job() );
			System.out.println("checking mutex now");
			if(mutex == 0)
			{	mutex = 1;}
			else
			{
				saveOptions(false);
				IJ.log("wrote " + options.getPath());
				IJ.error("There is already another active hexsampler Session!");
				return;
			}
			System.out.println(getClass().getSimpleName() + "\t*********NEW SESSION*********");
			IJ.showProgress(0.0);
			date.setTime(System.currentTimeMillis());
			IJ.log( "launching: " + dateFormat.format(date));
			if(monitor_worker >= 0)
			{
				System.out.println("forwarding stdout from: " + threadVec.get(monitor_worker).command);
			}
			write_results(false,true);
			saveOptions(false);
			worker_ready = true;
			boolean do_mainloop = true;
			new_submissions = true;
			int bt = 0;
			try
			{	
				do //mainloop
				{
					boolean new_state = false;
					while(bt < use_threads)
					{
						if( ( (busy_frames[bt]>=0) ) && (reader[bt].ready())) //dont even try to listen to zombie workers
						{	
							do_mainloop &= handleReply(bt);
							new_state = true;
						}
						if( (err_reader[bt].ready()))
						{	handleError(bt);} 
						++bt;
					}
					if(new_state)
					{	write_results(false,true);}
					bt %= use_threads;
					if(target_needs_update)
					{	
						update_uc_target();	
					}
					if( worker_ready /*|| ( new_submissions && (active_frames < use_frames) && (next_frame < use_frames)  )*/ )
					{
						IJ.showProgress(next_frame,use_frames);
						worker_ready = false;
						fill_busy_frames();
					}
					
					
					
					if(IJ.escapePressed())
					{
						IJ.resetEscape();
						IJ.showProgress(next_frame,use_frames);
						if(!IJ.showMessageWithCancel(getClass().getSimpleName()+ " halted",
														"active workers: " + active_frames +
														"\npending frames: " + (use_frames-next_frame) +
														"\ntotal frames done: " + frames_done +
														"\ntotal resets: " + total_resets +
														"\ndeaths/relaunches: " + deaths +
														"\ntotal delta merit: " + total_delta_merit +
														"\ntotal delta hex merit: " + total_delta_hex_merit +
													"\nClick Ok to Resume\nClick Cancel to " + 
													(new_submissions && endless?"release sisyphos" : 
													(sfImg==null?"flush samplers":"dumpsave subframes" )) + 
													"\nClick Cancel and hold Ctrl to Terminate") )
						{	
							if(IJ.controlKeyDown()) //the keystate is sometimes wrong
							{	
								IJ.log("Terminating hexsampler");
								do_mainloop = false;
							}
							
							if(new_submissions && !endless)
							{
								if(sfImg != null )
								{
									//collect_subframes(); //they should be up to date anyways
									dumpsave_sfImg();
								}
								else
								{
									if(new_submissions)
									{	new_submissions = false;}
									else
									{
										IJ.showStatus("flushing all workers");
										Process ps = rt.exec( IJ.getDirectory("imagej") + "flushsamplers.sh" );
										ps.waitFor();
										ps.destroy();
									}
									new_submissions = false;
								}
							}
							
							if(new_submissions && endless)
							{
								IJ.log("Sisyphos is sent to rest");
								new_submissions = false;
								forever = false; //TODO added this to see if that fixes shut down issue in endless mode
							}
							
							
						}	
					}
					synchronizeSlices(-1); //-1 determine slice from active imgage
					if( IJ.spaceBarDown() )
					{	showFrameInfo();	}
					
					if(do_mainloop) Thread.sleep(50);
				}
				while ( do_mainloop && (active_frames > 0) ); //mainloop
			}
			catch(Exception e)
			{
				IJ.log("Exception in mainloop");
				e.printStackTrace();
			}
			write_results(true,(op_mode!=2));//dont include filtered tasks when polling
			if(  (!batch_run &&  (op_mode == 3) ) || !do_mainloop )
			{
				if( IJ.showMessageWithCancel(getClass().getSimpleName() + " has finished",	
									"Clear empty subframes and sf_sums now?")
					)
				{	clear_empty_sf();}
			}
			else if (sfImg != null)
			{	clear_empty_sf();}
			
			shut_down();
			
			if(!batch_run && ( (debug_level > 0) || (op_mode == 3) ) )
			{	wait_for_Esc();} //user can inspect frames
			System.out.println(getClass().getSimpleName() + "\t*********END SESSION*********");
		} //for batch_id	
		mutex = 0; 
	}
	
	private void wait_for_Esc()
	{
		IJ.showMessage("Esc will quit " + getClass().getSimpleName() + "\n" +
					   "use space bar/ctrl and ROIs to inspect frames now\n");
		IJ.resetEscape();
		try
		{
			do
			{
				synchronizeSlices(-1); //-1 determine slice from current image
				if( IJ.spaceBarDown() )
				{	showFrameInfo();	}
				Thread.sleep(100);
			}	
			while( !IJ.escapePressed() );
		}
		catch(Exception e)
		{
			IJ.log("Exception in wait_for_Esc");
			e.printStackTrace();
		}
	}
	
	
	private void showFrameInfo()
	{
		ImagePlus imp = WindowManager.getCurrentImage();
		int slice = -1;
			
		if (imp != null)
		{
			if(imp == statImg || imp.getTitle().startsWith("Montage") || imp.getTitle().startsWith("Live_Montage"))
			{
				Roi roi = imp.getRoi();
				if ( roi != null )
				{
					Rectangle rect = roi.getBounds();
					int res = 1;
					if(imp.getTitle().startsWith("Montage") || imp.getTitle().startsWith("Live_Montage"))
					{
						
						res = uc_res; 
						String info = (String)imp.getProperty("Info");
						int mox = -1;
						try
						{	mox = Integer.parseInt( info.split("\n")[0].split("=")[1] );}
						catch (Exception e)
						{	e.printStackTrace();}
						if(mox > 0)
						{	res = imp.getWidth() / mox;}
						else
						{	IJ.error("Could not read xMontage of " + imp.getTitle() + "\n" + info);}
					}
					
					int x = rect.x/res;
					int y = rect.y/res;
					slice = 1 + x + y * (imp.getWidth()/res);
					//IJ.log("Reading Roi inside " + imp.getTitle() + " at " + rect.x + "," + rect.y + " as slice: " + slice);
				}
				else
				{	
					//IJ.log("Cursor lies outside " + imp.getTitle()); 
					return; //nothing to do
				}	
			}
			else 
			{	
				slice = imp.getCurrentSlice();
			}
			
			if (imp == Idata) //need to figure out which report (aka "slice") belongs to that frame
			{
				for( int i = 0; i < use_frames; ++i)
				{
					if( total_reports[i].frNr == slice)
					{
						slice = total_reports[i].repID;
						break;
					}		
				}
			}
			
			int rID = -1;//total_reports[slice-1].repID-1;//default for forward ordered frames
			for( int i = 0; i < use_frames; ++i)
			{
				if( total_reports[i].repID == slice)
				{
					rID = i;
					break;
				}		
				
			}
			if( (imp == statImg) || (imp == Idata) ||
			    (imp == uc_sum)  || (imp == uc_sym) ||
			    (imp == hexImg)  || (imp == sf_sum) ||
			    (imp == uc_mini) ||
			    imp.getTitle().startsWith("Montage") || imp.getTitle().startsWith("Live_Montage"))
			{
				if(rID >= 0 && rID < total_reports.length)
				{
						/*if(IJ.altKeyDown()) //alt is already in use by unity desktop
						{
							IJ.showStatus(	"Rp: " + rID +
											"Fr: "  + total_reports[rID].frNr	
											);
								
						}*/
						
					if(IJ.controlKeyDown())
					{
						IJ.showStatus(	"Rp: " + total_reports[rID].repID +
										" h_M: "  + total_reports[rID].hex_merit +
										" h_std: "  + total_reports[rID].hex_std +
										" contr(%): " + (100*total_reports[rID].contrast)	
										);	
					}
					else
					{
						IJ.showStatus("Fr: "  + total_reports[rID].frNr +
									 " M: "  + total_reports[rID].merit +
									 " runs: "   + total_reports[rID].runs +
									 " worker: " + total_reports[rID].worker +	
									 " " +  states[total_reports[rID].state+2] +
									 " t(s): " + total_reports[rID].elapsed);
					}
					if(imp == statImg || imp.getTitle().startsWith("Montage") || imp.getTitle().startsWith("Live_Montage") )
					{//synchronize the other stacks with the cursor position
						synchronizeSlices(total_reports[rID].repID);	
					}			 
					if(imp == Idata)
					{
						HexSampleReport report = total_reports[rID];
						double ox = report.offset_x + 0.5;
						double oy = report.offset_y + 0.5;
						
						
					
					
					
					}
				
				
				}
				
				else
				{	IJ.showStatus("Could not retrieve report " + rID);}
			}
			else if (imp == uc_target)
			{
				IJ.showStatus(" updates: " + target_updates + " in use: " + target_useable + " updating: " + update_target);
			}
		}	
	}
	
	private void synchronizeSlices( int set_slice)
	{
		boolean apply_slice = false;
		if(set_slice < 1) //retrieve actual slice from current image
		{
			ImagePlus imp = WindowManager.getCurrentImage();		
			if (imp != null)
			{
				int slice = imp.getCurrentSlice();
				//several reports may refer to the same Idata slice
				//so Idata is no key for synchronizing, but it is set if it is not virtual 
				
				if  ( 	(slice != the_slice) && 
						(	(imp == uc_mini) || (imp == Idata)   ||
							(imp == Mask)    || 
							(imp == uc_sum)  || (imp == uc_sym)  ||
							(imp == hexImg)  || (imp == sf_sum)
						)	
					) 
				{
					the_slice = slice;
					if(imp == Idata)
					{	
						for(int i = 0; i < total_reports.length; ++i)
						{
							//scan for matching frNR
							if(total_reports[i].frNr == the_slice) //frNr is one based
							{
								the_slice = total_reports[i].repID; //repID is also one based
								break;
							}
						}
					}
					else if (imp == Mask)
					{
						for(int i = 0; i < total_reports.length; ++i)
						{
							//scan for matching frNR
							if(total_reports[i].mkNr == the_slice) //mkNr is one based
							{
								the_slice = total_reports[i].repID; //repID is also one based
								break;
							}
						}
					
					}
					apply_slice = true;
				}
			}
		}
		else
		{
			the_slice = set_slice;
			apply_slice = true;
		}		
		if(apply_slice)
		{
			if(uc_sum != null)	uc_sum.setSlice(the_slice);
			if(uc_sym != null)	uc_sym.setSlice(the_slice);
			if(hexImg != null)	hexImg.setSlice(the_slice);
			if(sf_sum != null)	sf_sum.setSlice(the_slice);
			if(uc_mini != null)	uc_mini.setSlice(the_slice);
			if(Idata != null && !IdataSt.isVirtual())	
			{	
				for(int i = 0; i < total_reports.length; ++i)
				{
					//scan for matching repID
					if(total_reports[i].repID == the_slice) //repID is one based
					{
						Idata.setSlice(total_reports[i].frNr);
						break;
					}
				}
			}
			if(has_mask && !MaskSt.isVirtual())	
			{	
				for(int i = 0; i < total_reports.length; ++i)
				{
					//scan for matching repID
					if(total_reports[i].repID == the_slice) //repID is one based
					{
						Mask.setSlice(total_reports[i].mkNr);
						break;
					}
				}
			}
		}
	}
	
	private void clear_empty_sf()
	{
		int dead_sf_sums = 0;
		int dead_sfs = 0;
		
		if(sf_sumSt != null)
		{
			for(int sl = 1; sl <= sf_sumSt.getSize(); ++sl)
			{
				long sum = 0;
				short[] pixels = (short[])sf_sumSt.getPixels(sl);
				for(int i = 0; i < pixels.length; ++i)
				{	sum += (long)pixels[i];}
				if(sum == 0)
				{	
					if(sl < sf_sumSt.getSize())
					{	sf_sumSt.deleteSlice(sl);}
					else
					{	sf_sumSt.deleteLastSlice();}
					++dead_sf_sums;
					--sl;
				}
			}
			sf_sum.setStack(sf_sumSt);
			sf_sum.updateAndRepaintWindow();
			try
			{
				Thread.sleep(10);
			}
			catch( Exception e)
			{
				e.printStackTrace();
			}	
		}
		if( sfImgSt != null)
		{
			for(int sl = 1; sl <= sfImgSt.getSize(); ++sl)
			{
				long sum = 0;
				short[] pixels = (short[])sfImgSt.getPixels(sl);
				for(int i = 0; i < pixels.length; ++i)
				{	sum += (long)pixels[i];}
				if(sum == 0)
				{	
					if(sl < sfImgSt.getSize())
					{	sfImgSt.deleteSlice(sl);}
					else
					{	sfImgSt.deleteLastSlice();}
					++dead_sfs;
					--sl;	
				}
			}
			sfImg.setStack(sfImgSt);
			sfImg.updateAndRepaintWindow();
			try
			{
				Thread.sleep(10);
			}
			catch(Exception e)
			{
				e.printStackTrace();
			}	
		}
		if(dead_sf_sums > 0 || dead_sfs > 0)
		{
			IJ.log("deleted " + dead_sf_sums + " empty sf_sums and " + dead_sfs + " empty sfs");
		}
		else
		{
			IJ.log("There are neither empty sf_sums nor empty subframes");
		}
	}
	
	
	private boolean shut_down()
	{
		try
		{
			
			//saveOptions(false);
			for (int i=0; i<use_threads; ++i)
			{
				//if( busy_frames[i] != -1 )
				//{
					pwos[i].println("Run(0)");
					pwos[i].flush();	
				//}
			}
			if(streamnum > use_threads)
			{
				pwos[use_threads].println("Run(0)");
				pwos[use_threads].flush();
				pwos[use_threads].close();
				os[use_threads].close();
			}
			
			if(!keep_Idata) {Idata.close();}
			date.setTime(System.currentTimeMillis());
			IJ.log( "finished: " + dateFormat.format(date));
			IJ.log( "total delta merit: " + total_delta_merit + " total delta hex merit: " + total_delta_hex_merit);
			update_mean_parameters(-1);//-1 nobody is excluded
			IJ.log("global average tilt,hel,phi,excent: " +
					 mean_tilt + ", " + mean_hel + ", " + mean_phi + ", " + (mean_ellA/mean_ellB) );
			
			IJ.showProgress(1.0);
			/*for(int i=0; i < use_threads; ++i) 
			{	//This only floods the log
				IJ.log("Worker" + i + ": " + fr_p_worker[i] + " frames");
			}*/
			String prfx = prefix; 
			if(prefix.endsWith("_"))
			{	prfx = prefix.substring(0,prefix.length()-1);	}
			String logname = ("log_" + prfx + ".txt");
			
			IJ.log("saving " + logname);
			boolean log_ok = saveLog(path, logname);
			if(log_ok)
			{
				if(!batch_run)
				{	IJ.open(path + logname);}
				if(WindowManager.getWindow("Log") != null)
				{
					IJ.selectWindow("Log");
					IJ.run("Close");
				}
			}
			for (int i=0; i<use_threads; ++i)
			{	p[i].destroy(); }
			HexSampleReport.label_ok = true; //reenable the bad label warnings
			return true;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			IJ.log("Exception during shut_down");
			return false;
		}
	}
	
	private void handleError(int i)
	{
		try
		{
			boolean evil_error = false;
			if(err_reader[i].ready())
			{
				do
				{
					String line = err_reader[i].readLine();
					IJ.log( line );
					evil_error |= line.startsWith("EVIL");
				}
				while( err_reader[i].ready());			
				
				
				if(evil_error)
				{	
					//dont bother any longer
					threadVec.get(i).relaunches = max_relaunches; 
				}
				relaunch(i);
			}
		}
		catch( Exception e )
		{
			System.out.println("Exception while reading error stream of worker" + i);
			e.printStackTrace();
		}
	}
	
	private void relaunch(int i) 
	throws Exception
	{
		int rID = busy_frames[i];
		if(rID>=0)
		{
			IJ.log("dead worker" + i + " on Frame: " + total_reports[rID].repID );
			//pwos[i].println("Run(0)");
			//pwos[i].flush();
			//busy_frames[i]=-2;//disregard worker i in future
			--active_frames;
			total_reports[rID].worker = 0;
			total_reports[rID].state = 0; //pending
			p[i].destroy();
			++deaths;
			if(threadVec.get(i).relaunches < max_relaunches)
			{
				p[i] = rt.exec(threadVec.get(i).command);  //launch thread i
				os[i] = p[i].getOutputStream(); //get process i's input stream to os[i]
				pwos[i] = new java.io.PrintWriter(os[i]);
				is[i] = new BufferedInputStream( p[i].getInputStream() );  //get process i's output to is[i]
				dis[i] = new DataInputStream(is[i]);
				err_s[i] = p[i].getErrorStream();
				reader[i] = new java.io.InputStreamReader(dis[i]);
				err_reader[i] = new java.io.BufferedReader(new java.io.InputStreamReader(err_s[i]));
				busy_frames[i]=-1;
				worker_ready = true;
				threadVec.get(i).relaunches++;
				IJ.log("launched another process: " + threadVec.get(i).command);
			}
			else
			{
				IJ.log("giving up relaunching: " + threadVec.get(i).command);
			}
		}
	}
	
	private boolean handleReply(int i)
	{
		String line = "<null>";
		try
		{
			//boolean do_loop = true;
			/* //we no longer even listen to zombies
			if(busy_frames[i] == -1)
			{
				IJ.log("Illegal reply from zombie" + i + " dumping its stream now ...");
				while( reader[i].ready() )
				{	dis[i].readLine();}
				IJ.log("done");
				return true;
			}
			*/
			int rID = busy_frames[i];
			boolean final_report = false;
			boolean report_finished = false;
			boolean incomming_report = false;
			boolean first_warning = true;
			boolean new_uc_mini = false;
			double elapsed = 0.0;
			int report_lines = 100; 
			while( (--report_lines > 0) && ( incomming_report || reader[i].ready() ) )
			{
				if(report_lines == 1)
				{
					IJ.log("Error Report of worker" + i + " is too long");
				}
				
				line = dis[i].readLine();
				if(line.startsWith("REPORT"))
				{
					//IJ.log("requesting report from worker" + i);
					pwos[i].println("Report(1)");
					pwos[i].flush();
					final int wait_ms = 25; //that should be enough time to respond
					for(int w=0; w < 100; ++w) //max waiting time is 2*100*25ms = 5s
					{
						Thread.sleep(wait_ms);//most of the times 25ms are sufficient 
						if(reader[i].ready())
						{
							line = dis[i].readLine();
							incomming_report = true;
							break;
							//IJ.log("reading report of worker" + i);
						}
						Thread.sleep(wait_ms); 
					}
					if(!incomming_report)
					{
						IJ.log("worker"+i+" FAILED to deliever a promised report within 2.5s");
						relaunch(i);
						//technically the master is still fine and mainloop can shall go on
						return true;
					}
				}
				
				
				if ( ( (busy_frames[i] < 0) || (rID == -1) ) && !line.startsWith("Message"))
				{
					IJ.log("Zombie" + i + " moans: " + 	line);
				}
				if( (i == monitor_worker) && (debug_level > 1) ) //TODO check that
				{	System.out.println(line);}
				
				Scanner lsc = new Scanner(line);
				if(!lsc.hasNext())
				{	break;}
				String key = lsc.next();
				
				if( key.equals("Frame") )
				{
					int fr = lsc.nextInt();
					if( (fr != -1) && (fr != total_reports[rID].frNr) )
					{
						IJ.log("WARNING: Worker" + i + " reports on frame " + fr + " instead of " +  total_reports[rID].frNr);
						//busy_frames[i] = fr;
						total_reports[rID].frNr = fr;
						
					}
					if(fr == total_reports[rID].frNr)
					{	final_report = true;}	
				}
				else if( key.equals("HexEdgeLen") )
				{
					double hel = lsc.nextDouble();
					total_reports[rID].hexedgelen = hel;
				}
				else if( key.equals("Tilt") )
				{
					double tilt = lsc.nextDouble();
					total_reports[rID].tilt = tilt;
				}
				else if( key.equals("EllipseA") )
				{
					double ea = lsc.nextDouble();
					total_reports[rID].ellipseA = ea;
				}
				else if( key.equals("EllipseB") )
				{
					double eb = lsc.nextDouble();
					total_reports[rID].ellipseB = eb;
				}
				else if( key.equals("EllipsePhi") )
				{
					double phi = lsc.nextDouble();
					total_reports[rID].ellipsePhi = phi;
				}
				else if( key.equals("OffsetX") )
				{
					double dx = lsc.nextDouble();
					total_reports[rID].offset_x = dx;
				}
				else if( key.equals("OffsetY") )
				{
					double dy = lsc.nextDouble();
					total_reports[rID].offset_y = dy;
				}
				else if( key.equals("UC_Avg") )
				{
					double raw_avg = lsc.nextDouble();
					total_reports[rID].raw_avg = raw_avg;
				}
				else if( key.equals("Hex_Avg") )
				{
					double hex_avg = lsc.nextDouble();
					total_reports[rID].hex_avg = hex_avg;
				}
				else if( key.equals("Hex_Low") )
				{
					double hex_low = lsc.nextDouble();
					total_reports[rID].hex_low = hex_low;
				}
				else if( key.equals("Hex_High") )
				{
					double hex_high = lsc.nextDouble();
					total_reports[rID].hex_high = hex_high;
				}
				else if( key.equals("Hex_Merit") )
				{
					double hex_merit = lsc.nextDouble();
					total_reports[rID].hex_merit = hex_merit;
				}
				else if( key.equals("Hex_Shape") )
				{
					double hex_shape = Double.NaN; 
					if(lsc.hasNextDouble())
					{	hex_shape = lsc.nextDouble();}
					total_reports[rID].hex_shape = hex_shape;
				}
				else if( key.equals("Moment2") )
				{
					double moment2 = Double.NaN; 
					if(lsc.hasNextDouble())
					{	moment2 = lsc.nextDouble();}
					total_reports[rID].moment2 = moment2;
				}
				else if( key.equals("Contrast") )
				{
					double contrast = lsc.nextDouble();
					total_reports[rID].contrast = contrast;
				}
				else if( key.equals("Position_Q") )
				{
					double posQ = lsc.nextDouble();
					total_reports[rID].positionQ = posQ;
				}
				else if( key.equals("Mirror_Q") )
				{
					double mirQ = lsc.nextDouble();
					total_reports[rID].mirrorQ = mirQ;
				}
				else if( key.equals("SubFrames") )
				{
					int sf = lsc.nextInt();
					//check that there are as many subframes as expected
					//well we can handle any number anyways
					if(debug_level > 1)
					{
						if( (sfImg != null) && (sf != total_reports[rID].subframes) )
						{	IJ.log("Warning: changed number of subframes for frame " + total_reports[rID].frNr +
							 " reported: " + sf + " expected: " + total_reports[rID].subframes); }
					}
					total_reports[rID].subframes = sf;
				}
				else if( key.equals("OffsetHexQ") )
				{
					int dq = lsc.nextInt();
					total_reports[rID].offset_hexq = dq;
				}
				else if( key.equals("OffsetHexR") )
				{
					int dr = lsc.nextInt();
					total_reports[rID].offset_hexr = dr;
				}
				else if( key.equals("Merit") )
				{
					double merit = lsc.nextDouble();
					total_reports[rID].merit = merit;
				}
				else if( key.equals("sumUC") )
				{
					int area = lsc.nextInt();
					read_sum_uc(i, area);
				}
				else if( key.equals("symUC") )
				{
					int area = lsc.nextInt();
					read_sym_uc(i, area);
				}
				else if( key.equals("miniUC") )
				{
					int area = lsc.nextInt();
					new_uc_mini = read_mini_uc(i, area);
				}
				/*else if( key.equals("requestingTarget") )	
				{	
					//IJ.log("received target request from worker" + i);
					resubmit_target = true;
				}*/
				else if( key.equals("hexImg") )
				{
					int area = lsc.nextInt();
					read_hexImg(i, area);
				}
				else if( key.equals("sfSum") )
				{
					int area = lsc.nextInt();
					read_sf_sum(i, area);
				}
				else if( key.equals("sfStack") )
				{
					int area = lsc.nextInt();
					read_sf_stack(i, area);
				}
				else if( key.equals("Elapsed") )
				{
					elapsed = lsc.nextDouble();
					total_reports[rID].elapsed = elapsed;
					incomming_report = false;
					report_finished = true;		
				} 
				else if( key.equals("Message") )
				{
					String mess = lsc.nextLine();
					mess = mess.trim();
					if( (mess.startsWith("WARNING")) || (debug_level > 1) )
					{	IJ.log("worker" + i + ": " + mess);}
					//break;
				}
				else if(key.equals("EXCEPTION"))
				{
					incomming_report = false;
					total_reports[rID].state = 0;
					IJ.log("frame: " + total_reports[rID].frNr + " is withdrawn from worker" + i);
					--active_frames;
					busy_frames[i] = -1;
				}
				else if (i==monitor_worker)
				{
					if(debug_level < 2) //we already printed it at debug_level > 1
					{	System.out.println(line);} //break;
				}
				else if (first_warning)
				{
					IJ.log("Unexpected reply worker"+ i + ": " + line);
					IJ.log("Suppressing further errors from this read");
					first_warning = false;	
				}
				 //give the workers a chance to send all at once
			}
			if(report_finished)
			{
				//TODO run the test
				//MAYBE updating and using leeds to unstable feedback
				target_needs_update = (update_target && new_uc_mini); 
				if(final_report)
				{
					double delta = total_reports[rID].merit - total_reports[rID].old_merit;
					double hex_delta = total_reports[rID].hex_merit - total_reports[rID].old_hex_merit;
					total_delta_merit += delta;
					total_delta_hex_merit += hex_delta;
					/*
					if( (debug_level > 1) || 
						(total_reports[rID].hex_merit < 0) ||
						( ( (delta > 0.1*total_reports[rID].old_merit) || (hex_delta > 0.1*total_reports[rID].old_hex_merit) ) && 
						  ( (hex_delta>1) || (delta > 1) ) ) )
					{
						IJ.log("worker" + i + " frame: " + total_reports[rID].frNr + " hex_merit: " + total_reports[rID].hex_merit +
						" delta: " + hex_delta +
						" time: " + elapsed + "s" );
					}
					*/ 
					double rel_delta = hex_delta/total_reports[rID].hex_merit; //relative change may be negative in bad cases
					//For hex_merits < 10 we require at least 0.1 increase and 1% for larger values 
					total_reports[rID].state = ( forever && (hex_delta>0.1) && (rel_delta>0.01) ) ? 3 : -1; // 3.. finished and redo -1 finished and converged 
					if(forever && endless)
					{	
						total_reports[rID].state = 3;
					}
					
					//if(op_mode != 0) //true sisyphus mode
					//{	total_reports[rID].state = 3;}
					total_reports[rID].old_merit = total_reports[rID].merit;
					total_reports[rID].old_hex_merit = total_reports[rID].hex_merit;
					--active_frames;
					++frames_done;
					worker_ready = true;
					busy_frames[i] = -1;
					//hexedgelen should be quite uniform within a batch 10% deviaton triggers a reset 
					///Hmmm What shall we do if the very first frame failes ???
					if(	reset_bad_frames && forever && (op_mode==0) &&  ( total_reports[rID].merit < 30 ) && //better than 50 is typically super anyways
						update_mean_parameters(rID) &&
						( total_reports[rID].resets < 2) && 
						( 
						  //( total_reports[rID].contrast <= 0.025) || //min_hex_contrast) 						 
						  ( total_reports[rID].hex_merit < 0 ) || //this is definitely bad
						  ( total_reports[rID].hexedgelen < 0.96 * mean_hel) ||
						  (	total_reports[rID].hexedgelen > 1.04 * mean_hel) ||
						  (	total_reports[rID].offset_x < 0) ||
						  (	total_reports[rID].offset_y < 0) ||
						  (	total_reports[rID].offset_x > impWidth) ||
						  (	total_reports[rID].offset_y > impHeight) ||
						  ( Math.abs((total_reports[rID].tilt-mean_tilt)%(Math.PI/3)) > 0.05 ) ||
						  ( total_reports[rID].hex_shape < 1)  
					    )
					  )  
					{
						total_reports[rID].resets++;
						++total_resets;
						if(debug_level > 0)
						{	
							IJ.log("WARNING bad frame: " + total_reports[rID].frNr + " hex_merit: " + total_reports[rID].hex_merit + " tilt,hel,phi,excent: " +
							 total_reports[rID].tilt + ", " + total_reports[rID].hexedgelen + ", " + total_reports[rID].ellipsePhi + ", " + (total_reports[rID].ellipseA/total_reports[rID].ellipseB) );
							
							if(reset_bad_frames)
							{
								IJ.log("applying grouped N-1 average tilt,hel,phi,excent: " +
								mean_tilt + ", " + mean_hel + ", " + mean_phi + ", " + (mean_ellA/mean_ellB) );
							} 
							 
					    }
						
						if(reset_bad_frames)
						{
							apply_mean_parameters(rID);
							total_reports[rID].state = 0;
							total_reports[rID].runs = 0; //this will also reset the stability
							total_reports[rID].stability = 2;
						}
					}
					
					
				}
				else
				{
					if( total_reports[rID].state == 1) //fresh launch
					{
						if(total_reports[rID].runs == 1) //remember initial merit at first run
						{	
							total_reports[rID].old_merit = total_reports[rID].merit;
							total_reports[rID].old_hex_merit = total_reports[rID].hex_merit;	
						}
						total_reports[rID].state = 2; // now wait for final report
					}
					if(!IJ.spaceBarDown())
					{	IJ.showStatus("worker" + i + " frame: " +total_reports[rID].frNr +
						" hex_merit: " + total_reports[rID].hex_merit + "  status: " + states[total_reports[rID].state+2] );}
				}
				//This place is too busy
				//write_results(false,true);
				
			}
			return true;
		}
		catch(Exception e)
		{
			System.out.println("Exception while reading output stream of worker" + i + ": " + line);
			e.printStackTrace();
			return false;
		}	
	}
	
	private void read_sum_uc(int i, int len)
	//throws java.io.IOException
	{
		//IJ.log("worker" + i + " reading uc_sum");
		int repID = total_reports[busy_frames[i]].repID;
		short[] pixels = null;
		if( (uc_sum != null) && ( len == uc_sum.getWidth() *  uc_sum.getHeight()) )
		{ //read the data
			
			pixels = (short[])uc_sumSt.getPixels(repID);
		}
		else // do a fake read
		{
			IJ.log("commencing fake uc_sum read");
			pixels = new short[len];
		}
		if(!read_short_frame(i, pixels))
		{	
			System.out.println("error reading uc_sum" + repID + " from worker" + i);
			//throw new RuntimeException(Macro.MACRO_CANCELED);
		}
		if( uc_sum != null )	
		{	
			uc_sum.updateAndRepaintWindow();
			try
			{
				Thread.sleep(10);
			}
			catch(Exception e)
			{
				e.printStackTrace();
			}	
		}
	}
	
	private void read_sym_uc(int i, int len)
	//throws java.io.IOException
	{
		//IJ.log("worker" + i + " reading uc_pos");
		int repID = total_reports[busy_frames[i]].repID;
		short[] pixels = null;
		if( (uc_sym != null) && (len == uc_sym.getWidth() * uc_sym.getHeight() ) )
		{ //read the data
			pixels = (short[])uc_symSt.getPixels(repID);
		}
		else // do a fake read
		{
			IJ.log("commencing fake uc_sym read");
			pixels = new short[len];
		}
		if(!read_short_frame(i, pixels))
		{	
			System.out.println("error reading uc_sym" + repID + " from worker" + i);
			//throw new RuntimeException(Macro.MACRO_CANCELED);
		}
		if( uc_sym != null)
		{
			uc_sym.updateAndRepaintWindow();
			try
			{
				Thread.sleep(10);
			}
			catch(Exception e)
			{
				e.printStackTrace();
			}	
		}
	}
	
	private boolean read_mini_uc(int i, int len)
	//throws java.io.IOException
	{
		//IJ.log("worker" + i + " reading uc_pos");
		int repID = total_reports[busy_frames[i]].repID;
		short[] pixels = null;
		if( (uc_mini != null) && ( len == uc_mini.getWidth() * uc_mini.getHeight() ) )
		{ //read the data
			pixels = (short[])uc_miniSt.getPixels(repID);
		}
		else // do a fake read
		{
			IJ.log("commencing fake uc_mini read");
			pixels = new short[len];
		}
		if(!read_short_frame(i, pixels))
		{	
			System.out.println("error reading uc_mini" + repID + " from worker" + i);
			//throw new RuntimeException(Macro.MACRO_CANCELED);
			//return false;
		}
		if( uc_mini != null )	
		{
			uc_mini.updateAndRepaintWindow();
			try
			{
				Thread.sleep(10);
			}
			catch(Exception e)
			{
				e.printStackTrace();
			}	
			return true;
		}
		return false;
	}
	
	private void read_hexImg(int i, int len)
	//throws java.io.IOException
	{
		//IJ.log("worker" + i + " reading hex_Img");
		int repID = total_reports[busy_frames[i]].repID;
		short[] pixels = null;
		if( ( hexImg != null ) && ( len == hexImg.getWidth() * hexImg.getHeight() ) )
		{ //read the data
			pixels = (short[])hexImgSt.getPixels(repID);
		}
		else // do a fake read
		{
			IJ.log("commencing fake hexImg read");
			pixels = new short[len];
		}
		if(!read_short_frame(i, pixels))
		{	
			System.out.println("error reading hexImg" + repID + " from worker" + i);
			//throw new RuntimeException(Macro.MACRO_CANCELED);	
		}
		if( hexImg != null )
		{	
			hexImg.updateAndRepaintWindow();
			try
			{
				Thread.sleep(10);
			}
			catch(Exception e)
			{
				e.printStackTrace();
			}	
		}
	}
	
	private void read_sf_sum(int i, int len)
	//throws java.io.IOException
	{
		//IJ.log("worker" + i + " reading sf_sum");
		int repID = total_reports[busy_frames[i]].repID;
		short[] pixels = null;
		if((sf_sum != null) && ( len == sf_sum.getWidth() * sf_sum.getHeight() ) )
		{ //read the data
			pixels = (short[])sf_sumSt.getPixels(repID);
		}
		else // do a fake read
		{
			IJ.log("commencing fake sf_sum read");
			pixels = new short[len];
		}
		if(!read_short_frame(i, pixels))
		{	
			System.out.println("error reading sf_sum frame" + repID + " from worker" + i);
			//throw new RuntimeException(Macro.MACRO_CANCELED);	
		}
		if( sf_sum != null )
		{	
			sf_sum.updateAndRepaintWindow();
			try
			{
				Thread.sleep(10);
			}
			catch(Exception e)
			{
				e.printStackTrace();
			}	
		}
	}
	
	private void read_sf_stack(int i, int len)
	//throws java.io.IOException
	{
		//IJ.log("worker" + i + " reading sf");
		int rID = busy_frames[i];
		int repID = total_reports[rID].repID;
		String sl_lbl = sfImgArray[repID-1].getTitle();
		//int sfr = total_reports[busy_frames[i]].subframe_offset;
		int max_k = total_reports[rID].subframes;
		int last_k = max_k;
		boolean[] good_k = new boolean[max_k];
		short[][] pixels = new short[max_k][len];
		
		for(int k = 0; k < max_k; ++k)
		{
			good_k[k] = read_short_frame(i, pixels[k]);
			if(!good_k[k])
			{
				System.out.println("WARNING ignoring subframe (" + (k+1) + "of" + max_k + ") of frame " + repID + " from worker" + i);		
				--last_k;
				//throw new RuntimeException(Macro.MACRO_CANCELED);
			}	
		}
		total_reports[rID].subframes = last_k;
		ImageStack _sfImgSt = ImageStack.create(modelsize,modelsize,last_k,16);
		int kk = 0;
		for(int k = 0; k < max_k; ++k)
		{
			if(good_k[k])
			{
				++kk;
				_sfImgSt.setSliceLabel("(" + kk + "/" + last_k + ")" + sl_lbl, kk);
				_sfImgSt.setPixels(pixels[kk-1], kk);
			}
		}
		sfImgArray[repID-1].setStack(_sfImgSt);
		//sfImg.updateAndRepaintWindow();
	}
	
	private void collect_subframes()
	{
		int received_sf = 0;
		for(int fr = 0; fr < sfImgArray.length; ++fr)
		{
			final int sz = sfImgArray[fr].getImageStackSize();
			received_sf += (sz>1)?sz:0;
		}
		if(received_sf > 0)
		{
			sfImgSt = ImageStack.create(modelsize,modelsize,received_sf,16);
			int r_sf = 1;
			for(int fr = 0; fr < sfImgArray.length; ++fr)
			{
				ImageStack _sfImgSt = sfImgArray[fr].getStack();
				int sz = _sfImgSt.getSize();
				for(int sl = (sz>1)?1:2 ; sl <= sz; ++sl)
				{
					short[] pixels = (short[])_sfImgSt.getPixels(sl);
					String sl_lbl = _sfImgSt.getSliceLabel(sl);
					
					sfImgSt.setPixels(pixels, r_sf);
					sfImgSt.setSliceLabel(sl_lbl, r_sf);
						
					++r_sf;
				}
			}
			if(sfImg == null)//just in case the window was closed in the mean time
			{
				String d_title = ( !prefix.equals("") && Datatitle.startsWith(prefix) ) ? Datatitle.substring( prefix.length() ) : Datatitle;
				sfImg = NewImage.createShortImage( prefix + "sf_" + d_title, modelsize , modelsize , 1, NewImage.FILL_BLACK);
				sfImg.show();
			}
			sfImg.setStack(sfImgSt);
			sfImg.updateAndRepaintWindow();
			try
			{
				Thread.sleep(10);
			}
			catch(Exception e)
			{
				e.printStackTrace();
			}	
		}
	}
	
	
	
	
	private boolean read_short_frame(int i, short[] pixels)
	{
		try
		{
			//signal that we are ready to receive
			pwos[i].println("SendPixels(1)");
			pwos[i].flush();
			
			String line = dis[i].readLine();
			Scanner lsc = new Scanner(line);
			String key = lsc.next();
			if(key.equals("BeginPixels"))
			{
				final int len = lsc.nextInt();
				char[] buf = {0,0};
				int remaining = len;
				int skipped = 0;
				if(len == pixels.length)
				{
					long chck_sum = 0;
					for(int ind=0; ind < len; ++ind)
					{
						final int val = dis[i].readUnsignedShort();
						if (val == 65535)
						{
							pixels[ind] = 0;
							chck_sum += 0;
							System.out.println("Warning: Corrected MaxValue Pixel");
						} else 
						{			
							pixels[ind] = (short)val;
							chck_sum += val;
						}
					}
					line = dis[i].readLine();
					lsc = new Scanner(line);
					key = lsc.next();
					if(key.equals("EndPixels"))
					{
						long check = lsc.nextLong();
						if(check == chck_sum)
						{	
							if(check != 0)
							{	return true;}
							System.out.println("Error: received an empty image from worker " + i);
							//throw new RuntimeException(Macro.MACRO_CANCELED);
							return false;	
						}
						System.out.println("Error reading Pixels from "+ threadVec.get(i).command + ": " +
						 check + " != " + chck_sum + " : in " + getClass().getSimpleName() );
						//throw new RuntimeException(Macro.MACRO_CANCELED);
					}
					else
					{	
						System.out.println("Error reading from " + threadVec.get(i).command + ": " + key + " != EndPixels");
						//throw new RuntimeException(Macro.MACRO_CANCELED);	
					}
				}
				else
				{	
					System.out.println("Error reading from "+ threadVec.get(i).command + ": " + len + " != " + pixels.length );
					//throw new RuntimeException(Macro.MACRO_CANCELED);	
				}
			}
			else
			{	
				System.out.println("Error reading from "+ threadVec.get(i).command + ": " + line + " != BeginPixels " + pixels.length );
				//throw new RuntimeException(Macro.MACRO_CANCELED);	
			}
			return true;
		}
		catch( Exception e)
		{
			e.printStackTrace();
			return false;
		}
	}

	private boolean write_results(boolean last, boolean all_reports)
	{
		try
		{
			String prfx = prefix; 
			if(prefix.endsWith("_"))
			{	prfx = prefix.substring(0,prefix.length()-1);	}
			String[] fr_names = {"frame_continue.txt", "frame_" + prfx + ".txt"};
			
			for(int i = 0; (i < fr_names.length - (last?0:1)); ++i )
			{
				java.io.FileOutputStream dumpstream = new java.io.FileOutputStream(path + fr_names[i]);
				java.io.PrintWriter dumpwriter = new java.io.PrintWriter(dumpstream,false);
				date.setTime(System.currentTimeMillis());
				dumpwriter.println("#" + dateFormat.format(date) );
				dumpwriter.println("#" + Datatitle + "   " + use_frames);
				
				int last_match = -1;
				int last_group = -1;
				for(int rID=1;rID<=total_reports.length;++rID)
				{
					int searches = 0;
					do
					{
						int first_ff = 0;
						int last_ff = total_reports.length;
						if(last_match != -1 )
						{
							if(last_match > first_ff + 1)	first_ff = last_match - 1;
							if(last_match < last_ff-1) last_ff = last_match + 1;
							last_match = -1;
						}
						for(int ff=first_ff;ff<last_ff;++ff)
						{
							if( total_reports[ff].repID == rID)	
							{
								if( (last_group != total_reports[ff].group) && (total_reports[ff].group != -1) ) //groups must be actually grouped
								{
									last_group = total_reports[ff].group;
									dumpwriter.println("");
									if( last_group < groupHeaders.size() )
									{	dumpwriter.print(groupHeaders.get(last_group)); }
									dumpwriter.println(HexSampleReport.headerline);
									dumpwriter.println("");
								}
								last_match = ff;
								
								
								if(all_reports || (total_reports[ff].state != -2) )//skipped frames would mess up the counting of subframes
								{	dumpwriter.println(total_reports[ff].toString());}
								
								break;
							}
						}
					}
					while( (last_match == -1) && (++searches < 2) );
				}
				dumpwriter.close();
			}
			if(last)
			{	IJ.log("wrote " + path + fr_names[fr_names.length-1]);}
			
			if(sfImg != null)
			{	
				collect_subframes();
				if( (max_stack_size > 0) && (sfImg.getStackSize() >= max_stack_size))
				{	dumpsave_sfImg();}
			}
			
			if( (statImg != null) && (statImgSt != null) )
			{
				float[] merits = (float[])statImgSt.getPixels(1);
				float[] numSfs = (float[])statImgSt.getPixels(2);
				float[] raw_avgs = (float[])statImgSt.getPixels(3);
				float[] hex_avgs = (float[])statImgSt.getPixels(4);
				float[] hex_lows = (float[])statImgSt.getPixels(5);
				float[] contrasts = (float[])statImgSt.getPixels(6);
				float[] raw_stds = (float[])statImgSt.getPixels(7);
				float[] hex_stds = (float[])statImgSt.getPixels(8);
				float[] tilts = (float[])statImgSt.getPixels(9);
				float[] hels = (float[])statImgSt.getPixels(10);
				float[] phis = (float[])statImgSt.getPixels(11);
				float[] excents = (float[])statImgSt.getPixels(12);
				float[] positionQs = (float[])statImgSt.getPixels(13);
				float[] mirrorQs = (float[])statImgSt.getPixels(14);
				float[] shapeQs = (float[])statImgSt.getPixels(15);
				float[] timings = (float[])statImgSt.getPixels(16);
				float[] hex_highs = (float[])statImgSt.getPixels(17);
				float[] hex_merits = (float[])statImgSt.getPixels(18);
				float[] hex_shapes = (float[])statImgSt.getPixels(19);
				float[] syms = (float[])statImgSt.getPixels(20);
				float[] moments = (float[])statImgSt.getPixels(21);
				
				for(int ff=0; ff<total_reports.length; ++ff)
				{
					//int fr = total_reports[ff].frNr - 1;
					int repID = total_reports[ff].repID-1;
					merits[repID] = (float)total_reports[ff].merit;
					numSfs[repID] = (float)total_reports[ff].subframes;
					raw_avgs[repID] = (float)total_reports[ff].raw_avg;
					hex_avgs[repID] = (float)total_reports[ff].hex_avg;
					hex_lows[repID] = (float)total_reports[ff].hex_low;
					contrasts[repID] = (float)total_reports[ff].contrast;
					raw_stds[repID] = (float)( total_reports[ff].raw_std ); 
					hex_stds[repID] = (float)( total_reports[ff].hex_std );
					tilts[repID] = (float)total_reports[ff].tilt;
					hels[repID] = (float)total_reports[ff].hexedgelen;
					phis[repID] = (float)total_reports[ff].ellipsePhi;
					excents[repID] = (float)(total_reports[ff].excent);
					positionQs[repID]=(float)total_reports[ff].positionQ;
					mirrorQs[repID]=(float)total_reports[ff].mirrorQ;
					shapeQs[repID]=(float)total_reports[ff].shapeQ;
					timings[repID]=(float)total_reports[ff].elapsed;
					hex_highs[repID]=(float)total_reports[ff].hex_high;
					hex_merits[repID]=(float)total_reports[ff].hex_merit;
					hex_shapes[repID]=(float)total_reports[ff].hex_shape;
					syms[repID]=(float)total_reports[ff].sym;
					moments[repID]=(float)total_reports[ff].moment2;
				}	
				statImg.updateAndRepaintWindow();
				Thread.sleep(10);
			}
			return true;
		}
		catch(Exception e)
		{
			//resultwriter.flush();
			e.printStackTrace();
			return false;
		}	
		
	}
	
	
	private boolean initialize_job()
	{
		try 
		{
			if(!employThreads(threadnum))
			{	
				releaseThreads(use_threads);
				return false;
			}
			if(!init_frame_stats()) //initalizes from frame_file 
			{	
				releaseThreads(use_threads);
				return false;
			}
			saveOptions(false); ///Sometimes the JVM freezes in the code below !!!
			p = new java.lang.Process[max_threads];
			//the process' output (our input stream- here we read responses from the worker thread)
			is = new java.io.BufferedInputStream[max_threads];
			dis = new java.io.DataInputStream[max_threads];
			reader = new java.io.InputStreamReader[max_threads];
			err_s = new java.io.InputStream[max_threads];
			err_reader = new java.io.BufferedReader[max_threads];
			//the process' input stream (our output - here we send data and commands to it)
			os = new java.io.OutputStream[max_threads + 1];
			pwos = new java.io.PrintWriter[max_threads + 1];
			
			for (int i = 0; i < use_threads; ++i)
			{
				//int worker = threadVec.get(i).worker + thread_offset;
				p[i] = rt.exec(threadVec.get(i).command);  //launch thread i
				os[i] = p[i].getOutputStream(); //get process i's input stream to os[i]
				pwos[i] = new java.io.PrintWriter(os[i]);
				is[i] = new BufferedInputStream( p[i].getInputStream() );  //get process i's output to is[i]
				dis[i] = new DataInputStream(is[i]);
				reader[i] = new java.io.InputStreamReader(dis[i]);
				err_s[i] = p[i].getErrorStream();
				err_reader[i] = new java.io.BufferedReader(new java.io.InputStreamReader(err_s[i]));	
			}
			streamnum = use_threads;
			if(debug_level > 1)
			{
				++streamnum;
				IJ.log("writing " + path + prefix + "job.bin");
				java.io.FileOutputStream jobstream = new java.io.FileOutputStream(path + prefix + "job.bin", false);
				os[use_threads] =  new java.io.DataOutputStream(jobstream);
				pwos[use_threads] = new java.io.PrintWriter(jobstream);
			}
			/*
			if(uc_target_title.equals("<new>"))
			{
				uc_target = NewImage.createShortImage( prefix + "uc_target" , 2*minilength , 2*minilength , 1, NewImage.FILL_BLACK);
				uc_target_title = uc_target.getTitle();
				uc_target.show();
				target_useable = false;
				update_target = true;
			}
			
			if(uc_target != null)
			{	uc_targetSt = uc_target.getStack();}
			*/
			
			if( (op_mode == 0)  || (op_mode == 3) || (op_mode == 4))
			{
				uc_sum = NewImage.createShortImage( prefix + "uc_avg" , uc_res , uc_res , use_frames, NewImage.FILL_BLACK);
				uc_sumSt = uc_sum.getStack();
				uc_sym = NewImage.createShortImage( prefix + "uc_pos" , uc_res , uc_res , use_frames, NewImage.FILL_BLACK);
				uc_symSt = uc_sym.getStack();
				uc_mini = NewImage.createShortImage( prefix + "uc_mini" , 2*minilength , 2*minilength , use_frames, NewImage.FILL_BLACK);
				uc_miniSt = uc_mini.getStack();
				for(int f = 0; f < use_frames; ++f)
				{	
					String lbl = IdataSt.getSliceLabel(total_reports[f].frNr);
					int i = total_reports[f].repID;
					uc_sumSt.setSliceLabel(lbl,i);
					uc_symSt.setSliceLabel(lbl,i);
					uc_miniSt.setSliceLabel(lbl,i);
				}
				uc_sum.setStack(uc_sumSt);
				uc_sym.setStack(uc_symSt);
				uc_mini.setStack(uc_miniSt);
				Thread.sleep(10);
				uc_sum.show();
				Thread.sleep(10);
				uc_sym.show();
				Thread.sleep(10);
				uc_mini.show();
				Thread.sleep(10);
			}
			
			//should also be fine without the first check against "" 
			String d_title = ( !prefix.equals("") && Datatitle.startsWith(prefix) ) ? Datatitle.substring( prefix.length() ) : Datatitle;
			
			
			if( (op_mode == 1) || (debug_level > 0) )
			{
				int hex_width = (int)(impWidth*hexImgScale);
				int hex_height = (int)(impHeight*hexImgScale);
				if( hex_width % 2 == 1 ) {++hex_width;}
				if( hex_height % 2 == 1 ) {++hex_height;}	
				
				hexImg = NewImage.createShortImage( prefix + "hexed_" + d_title , hex_width , hex_height , use_frames, NewImage.FILL_BLACK);
				hexImgSt = hexImg.getStack();
				sf_sum = NewImage.createShortImage( prefix + "sf_sum", modelsize , modelsize , use_frames, NewImage.FILL_BLACK);
				sf_sumSt = sf_sum.getStack();
				for(int f = 0; f < use_frames; ++f)
				{	
					String lbl = IdataSt.getSliceLabel(total_reports[f].frNr);
					int i = total_reports[f].repID;
					hexImgSt.setSliceLabel(lbl,i);
					sf_sumSt.setSliceLabel(lbl,i);	
				}
				hexImg.setStack(hexImgSt);
				sf_sum.setStack(sf_sumSt);
				Thread.sleep(10);
				hexImg.show();
				Thread.sleep(10);
				
				sf_sum.show();
				//System.out.println("showed sf_sum"); //displayed if IJ crashed
				Thread.sleep(10);
			}
			
			
			if( (op_mode == 2) || ( (op_mode == 3) && (debug_level>=1) ) )
			{
				sfImgArray = new ImagePlus[use_frames];
				for(int f = 0; f < use_frames; ++f)
				{
					int sl = total_reports[f].frNr;
					String sl_lbl = IdataSt.getSliceLabel(sl);
					//There will be duplicate labels if there are several tasks per frame
					//will produce correct Warnings in Hex_Collector.java
					sfImgArray[f] = NewImage.createShortImage( prefix + sl_lbl, modelsize , modelsize , 1, NewImage.FILL_BLACK);
				}
				sfImg = NewImage.createShortImage( prefix + "sf_" + d_title, modelsize , modelsize , 1, NewImage.FILL_BLACK);
				sfImgSt = sfImg.getStack();
				sfImg.show();
				Thread.sleep(10);
				
				
			}
			
			{   //allways good to have the stat image
				
				int columns = (int)Math.sqrt(use_frames); //maybe we have to adjust this once we get much more than 400 frames  
				if(columns < 5) {	columns = 4; } 
				int height = ( use_frames + columns - 1 )/columns; 
				statImg = NewImage.createFloatImage( prefix + "stats_" + d_title, columns , height , 21, NewImage.FILL_BLACK);
				statImgSt = statImg.getStack();
				statImgSt.setSliceLabel("Merit",1);
				statImgSt.setSliceLabel("NumSF",2);
				statImgSt.setSliceLabel("Avg(uc)",3);
				statImgSt.setSliceLabel("Avg(hex)",4);
				statImgSt.setSliceLabel("Hex_low",5);
				statImgSt.setSliceLabel("Contrast(std/mean)",6);
				statImgSt.setSliceLabel("Std(uc)",7);
				statImgSt.setSliceLabel("Std(hex)",8);
				statImgSt.setSliceLabel("Tilt",9);
				statImgSt.setSliceLabel("HEL",10);
				statImgSt.setSliceLabel("Phi",11);
				statImgSt.setSliceLabel("Excent",12);
				statImgSt.setSliceLabel("Localization",13);
				statImgSt.setSliceLabel("MirrorQuality",14);
				statImgSt.setSliceLabel("ShapeQuality",15);
				statImgSt.setSliceLabel("Elapsed",16);
				statImgSt.setSliceLabel("Hex_high",17);
				statImgSt.setSliceLabel("Hex_merit",18);
				statImgSt.setSliceLabel("Hex_shape",19);
				statImgSt.setSliceLabel("Sym",20);
				statImgSt.setSliceLabel("Hex_Mom2",21);
				for(int sl = 1; sl <= statImgSt.getSize(); ++sl )
				{
					float[] pixels = (float[])statImgSt.getPixels(sl);
					for(int ind = 0; ind < pixels.length; ++ind)
					{	pixels[ind] = Float.NaN;}
				}
				statImg.setStack(statImgSt);
				statImg.show();
				Thread.sleep(10);
			}
			
			System.out.println("jobs were initialized");
			return true;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			IJ.error("Error during initialization. Please check for hexsampler zombies");
			releaseThreads(use_threads);
			return false;
		}
	}
	
	
	private boolean init_frame_stats()
	{
		try
		{
			boolean first_warning = true;
			boolean reply = true;
			boolean rerun = false;
			busy_frames = new int[use_threads];
			for(int i=0; i < use_threads; ++i)
			{	busy_frames[i]=-1; } //-1 .. no assigned report/task
			//busy_reports = new int[use_threads];
			fr_p_worker = new int[use_threads];
			
			do
			{
				rerun = false;
				ImagePlus IdataO = null;
				ImageStack IdataStO = null;
				ImagePlus MaskO = null;
				ImageStack MaskStO = null;
				
				
				total_reports = new HexSampleReport[use_frames];
				int avg_counter = 0;
				boolean ffok = false;
				Scanner ffs = null;
				groupHeaders.clear();
				//groupHeaders.add(new StringBuilder(256));
				int headergroup = -1;
				boolean header_core = false; //Wheter or not the original header created by Anreass ws encountered
				if(frame_file != null)
				{	
					ffs = new Scanner(new FileReader(frame_file));
					IJ.log("reading " + frame_file.getPath());
					ffok = (frame_file.isFile() && frame_file.setReadable(true)); 
					//should always be true unless file was deleted in the meantime
				}
			
				if(!ffok)
				{	
					IJ.log(" applying defaults hel: " + hexedgelen + "  tilt: " + tilt + " phi: " + def_phi + " excent(a/b): " + def_excent);
					frame_file = null;	
				}
				else if (reorder_stack)
				{
					IdataO = NewImage.createShortImage( prefix + "ordered_" + Datatitle , impWidth , impHeight , use_frames, NewImage.FILL_BLACK);
					IdataStO = IdataO.getStack();
					if(has_mask)
					{
						MaskO = NewImage.createShortImage( prefix + "ordered_" + Masktitle , impWidth , impHeight , use_frames, NewImage.FILL_BLACK);
						MaskStO = MaskO.getStack();
					}
				}
				HexSampleReport.def_phi = def_phi;
				HexSampleReport.def_excent = def_excent;
				HexSampleReport.current_group = -1; //groups of entries are incremented by set_key_lines
				if(filter_file != null)
				{	filter = new FrameFilter(filter_file);}
				int filtered = 0;
				final double rr = 2*impWidth/(3*hexedgelen*bondlength*Math.sqrt(3.0));
				for(int f = 0; f < use_frames; ++f)
				{	
					//framestates[i] = FrameState.PENDING;
					total_reports[f] = new HexSampleReport(hexedgelen, rr,tilt,f+1);
					if(ffok && ffs.hasNextLine())
					{
						String line = ffs.nextLine();
						//String header = line;
						line = line.trim();
						while (ffs.hasNextLine() && (line.startsWith("#") || line.length() == 0) )
						{
							if( line.startsWith("#label") || line.startsWith("#filename") )
							{	
								//header = line;
								HexSampleReport.set_key_line(line.substring(1));
								//System.out.println(line.substring(1));	
							}
							else if ( line.length() > 0 && (header_core || line.startsWith("#Information") || line.startsWith("#This") ))
							{
								header_core = true;
								while(groupHeaders.size() < headergroup+2)//typically a one time loop unless there are missing headers
								{
									groupHeaders.addElement( new StringBuilder(256) );
									//System.out.println("created new HeaderGroup <" + (headergroup + 1) + ">");	
								}
								groupHeaders.get(headergroup+1).append(line + "\n");
							}
							line = ffs.nextLine();
							line = line.trim();	
						}
						if( ! total_reports[f].parseString(line, Idata, IdataStO, MaskStO, reorder_stack, Mask) )
						{	return false;}
						/*
						if(headergroup != total_reports[f].group)
						{
							System.out.println("Header <" + (headergroup+1) + "> current Group: " + total_reports[f].group);
							System.out.println(groupHeaders.get(headergroup+1));
						}
						*/
						header_core = false; //for sure the headers must have finished AFTER a report was parsed
						headergroup = total_reports[f].group; //initially zero and incremented with every new #label
						total_reports[f].use_mirror = use_mirrorQ;
						total_reports[f].elapsed = 0.0; //reset the timers;
						total_reports[f].stability = stability;
						if(reset_states)
						{	total_reports[f].state = 0;}
						if( (filter != null) && !filter.apply_filter(total_reports[f]) )
						{
							++filtered;
							total_reports[f].state = -2; //skip the filtered reports
						}
						
						
					}
					else if (ffok)
					{	
						if(filtered == 0)
						{	IJ.log("Warning missing entry " + (f+1) + " (and beyond) in " + frame_file.getPath());}
						total_reports[f].state = -2; //skip the filtered reports
						++filtered; 	
					} 
				}
				if(filtered != 0)
				{	IJ.log("removed/filtered " + filtered + " out of " + use_frames + " frames(sampling tasks)");}
				if(ffok) 
				{	
					ffs.close();
					if(reorder_stack)
					{	
						if(has_mask)
						{
							IdataO.setStack(IdataStO);
							MaskO.setStack(MaskStO);
						}
						//IdataO.show();
						//IdataO.updateAndRepaintWindow();
						///We assume no idiot would ever change the order
						///of entries in frame.txt when combining different filters
						reorder_stack = false;
						
						Idata.close();
						Thread.sleep(10);
						Idata = IdataO;
						Datatitle = Idata.getTitle();
						Idata.show();
						Thread.sleep(10);
						IdataSt = Idata.getStack();
						impDepth = IdataSt.getSize();
						
						if(has_mask)
						{
							Mask.close();
							Thread.sleep(10);
							Mask = MaskO;
							Masktitle = Mask.getTitle();
							Mask.show();
							Thread.sleep(10);
							MaskSt = Mask.getStack();
							mskDepth = MaskSt.getSize();
						}
					}
					if(filtered != 0)
					{
						if( IJ.showMessageWithCancel(getClass().getSimpleName() + " asks"," Would you like to delete the filtered frames\n" + 
						 "and rewrite/reload frame_continue\n"+
						 "and reload " + filter_name + " now?\n" +
						 "You may modify " + filter_name + " before proceeding.") )
						{
							int[] used_frames = new int[use_frames];
							int used = 0;
							for(int f = 0; f < use_frames; ++f)
							{
								if(total_reports[f].state != -2)
								{
									if(used_frames[total_reports[f].frNr-1] == 0)
									{	++used;	} //several reports may refer to the same slice
									++used_frames[total_reports[f].frNr-1];
								}
								
							}
							IJ.log("There are " + used + " distinct frames in " + Idata.getTitle());
							
							ImageStack IdataStR = ImageStack.create(Idata.getWidth(),Idata.getHeight(),used,16);
							ImageStack MaskStR = null;
							if(has_mask)
							{
								MaskStR = ImageStack.create(Mask.getWidth(),Mask.getHeight(),used,8);
							}
							int g = 0;
							for(int f = 0; f < use_frames; ++f)
							{
								if(used_frames[f] != 0)
								{
									++g;//slices start from 1
									IdataStR.setPixels( IdataSt.getPixels(f+1), g);
									IdataStR.setSliceLabel(	IdataSt.getSliceLabel(f+1), g);
									if(has_mask)
									{
										MaskStR.setPixels( MaskSt.getPixels(f+1), g);
										MaskStR.setSliceLabel(	MaskSt.getSliceLabel(f+1), g);
									}
								}
							}
							Idata.setStack(IdataStR);
							IdataSt = Idata.getStack();
							Idata.updateAndRepaintWindow();
							Thread.sleep(10);
							impDepth = used;
							if(has_mask)
							{
								Mask.setStack(MaskStR);
								MaskSt = Mask.getStack();
								Mask.updateAndRepaintWindow();
								Thread.sleep(10);
								mskDepth = used;
							}
							write_results(false,false);//only write not skipped results to only frame_continue.txt
							//dont change the actual inputname, only the frame_file
							frame_file = new File(path + "frame_continue.txt");
							use_frames -= filtered; //number of filtered reports
							IJ.log("wrote " + (use_frames) + " tasks to " + frame_file.getPath());
							total_reports = null; //they are outdated
							rerun = true; 
						}
					}
						
				}
			}while(rerun);	
			if(backwards)
			{
				for(int f = 0; f < total_reports.length / 2; f++)
				{
					HexSampleReport temp = total_reports[f];
					total_reports[f] = total_reports[total_reports.length - f - 1];
					total_reports[total_reports.length - f - 1] = temp;
				}
			}
			
			update_mean_parameters(-1);//-1 nobody is excluded
			IJ.log("global average tilt,hel,phi,excent: " +
					 mean_tilt + ", " + mean_hel + ", " + mean_phi + ", " + (mean_ellA/mean_ellB) );
			System.out.println("Frames/Tasks were prepared");
			return reply;	
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return false;
		}	
	}
	
	private boolean update_mean_parameters(int i)
	{
		double new_mean_tilt_cos = 0.0;
		double new_mean_tilt_sin = 0.0;
		double new_mean_hel = 0.0;
		double new_mean_phi_cos = 0.0;
		double new_mean_phi_sin = 0.0;
		double new_mean_ellA = 0.0;
		double new_mean_ellB = 0.0;
		int num_pts = 0;
		double weight = 0.0;
		int home_group = -1;
		if(i >= 0 && i < use_frames)
		{	home_group = total_reports[i].group;}
		for(int f = 0; f < use_frames; ++f)
		{
			double hm = total_reports[f].hex_merit;
			double ctr = total_reports[f].contrast;
			double hs = total_reports[f].hex_shape; 
			double hm_2 = Math.pow(hm,2);
			if( (f != i) && ( (i<0) || (i >= use_frames) || (total_reports[f].group == home_group) ) &&
			(total_reports[f].state != -2) && (hs > 0.0) && (hm > 0.0) && (ctr > min_hex_contrast) )
			{
				double ellPhi = total_reports[f].ellipsePhi;
				double tilt = total_reports[f].tilt;
				ellPhi = (ellPhi%(2*Math.PI) + 2*Math.PI) % Math.PI;
				//tilt = (tilt%(2*Math.PI)+2*Math.PI) % (Math.PI/3);
				tilt = tilt % (Math.PI/3);
				new_mean_tilt_cos += hm_2 * Math.cos(6*tilt);
				new_mean_tilt_sin += hm_2 * Math.sin(6*tilt);
				new_mean_hel += hm_2 * total_reports[f].hexedgelen; //redundant
				new_mean_phi_cos += hm_2 * Math.cos(2*ellPhi);
				new_mean_phi_sin += hm_2 * Math.sin(2*ellPhi); 
				new_mean_ellA += hm_2 * total_reports[f].ellipseA;
				new_mean_ellB += hm_2 * total_reports[f].ellipseB;
				weight += hm_2;
				++num_pts;
			}
		}
		if(num_pts > 0)
		{
			double norm = 1.0/weight;
			//Math.atan2(y,x) gives phi
			mean_tilt = Math.atan2(new_mean_tilt_sin * norm, new_mean_tilt_cos * norm)/6;
			if(mean_tilt  > Math.PI/6)
			{	mean_tilt -= Math.PI/3;}
			if(mean_tilt <= -Math.PI/6)
			{	mean_tilt += Math.PI/3;}
			
			mean_hel = new_mean_hel * norm; //redundant
			mean_phi = Math.atan2(new_mean_phi_sin * norm, new_mean_phi_cos * norm)/2;
			if(mean_phi > Math.PI/2)
			{	mean_phi -= Math.PI;}
			if(mean_phi <= -Math.PI/2)
			{	mean_phi += Math.PI;}
			
			mean_ellA = new_mean_ellA * norm;
			mean_ellB = new_mean_ellB * norm;
			return true;
		}
		return false;
	}
	
	private boolean apply_mean_parameters(int i)
	{
		if( (i >= 0) && (i < use_frames) )
		{
			total_reports[i].tilt =	mean_tilt;
			total_reports[i].hexedgelen = mean_hel;	//redundant	
			total_reports[i].ellipsePhi = mean_phi;
			total_reports[i].ellipseA = mean_ellA;
			total_reports[i].ellipseB = mean_ellB;
			total_reports[i].excent = mean_ellA / mean_ellB;
			//total_reports[i].use_mirror = !total_reports[i].use_mirror;
			return true;
		}
		return false;	
	}
	
	
	private boolean fill_busy_frames()
	{
		boolean rerun = true;
		boolean shut_down_worker = false;
		for(int i = 0; i < use_threads; ++i)
		{
			if(busy_frames[i]==-1) //We found a free slot
			{
				next_frame = 0; //search the whole range
				while( (next_frame < use_frames) && (total_reports[next_frame].state != 0) && new_submissions ) // 0 .. pending
				{	++next_frame;}
				if(next_frame >= use_frames) //no more pending frames 
				{
					if( forever && new_submissions )
					{
						///The requeeuing seems not to kill workers
						int requeued = 0;
						for(int j = 0; j < use_frames; ++j)
						{	
							if(total_reports[j].state == 3) //done and ready for redo
							{	
								total_reports[j].state = 0; 
								++requeued;
								
							}//now pending
							//lets see if that will help
							if( (total_reports[j].resets > 0) && (frames_done < 2 * use_frames) )
							{
								total_reports[j].resets = 0;
								if(total_reports[j].state != 0) 
								{
									++requeued;
									total_reports[j].state = 0;
								}	
							} 	
							
						}
						if(!endless)
						{
							IJ.log("Sisyphos requeued frames: " + requeued);
						}
						next_frame = 0; //search the whole range
						while( (next_frame < use_frames) && (total_reports[next_frame].state != 0) && new_submissions) // 0 .. pending
						{	++next_frame;}
						if(next_frame >= use_frames) //even no more reactivated frames 
						{	shut_down_worker = true;}
					}
					
					if( !(forever && new_submissions) || shut_down_worker) 
					{
						shut_down_worker = false;
						if( i == monitor_worker)
						{	System.out.println("shutting down worker" + i);}
						pwos[i].println("Run(0)");
						pwos[i].flush();
						busy_frames[i]=-2;
						forever = endless;
					}	
				}
				
				if ( (next_frame < use_frames) && (total_reports[next_frame].state == 0) && new_submissions ) //recheck after changing the state
				{
					///DEBUG lets see if this suppresses deaths
					/// Well it seems to suppress the probability with 20
					/// lets observe with 25 for some time
					try
					{
						Thread.sleep(25);//maybe this is the place where workers get killed by "overfeeding"
					}
					catch (Exception e)
					{
						//no serious problem simply report and continue
						e.printStackTrace();
					}
					
					busy_frames[i] = (next_frame);
					//busy_reports[i] = next_frame;
					total_reports[next_frame].state = 1; // 1..busy
					total_reports[next_frame].old_merit = total_reports[next_frame].merit;
					total_reports[next_frame].old_hex_merit = total_reports[next_frame].hex_merit;
					boolean task_ok = submit_task(i);
					if((i==monitor_worker) && (streamnum > use_threads))
					{	
						if(!submit_task(use_threads))
						{	IJ.log("Error while writing to job.bin");}
					}
					if(task_ok)
					{	
						total_reports[next_frame].worker = i;
						++total_reports[next_frame].runs;	
					}
					else
					{	IJ.log("Error wile submitting frame " + total_reports[next_frame].frNr + " to worker" + i);}
						
					
					if(!task_ok) 
					{
						try
						{	
							if(err_reader[i].ready())
							{	IJ.log( err_reader[i].readLine() );	}
							relaunch(i);
						}
						catch(Exception e)
						{
							e.printStackTrace();
							return false;
						}
					}
				}
				
			}
		}
		return true;
	}
	

	private boolean submit_task(int i)
	{
		try
		{
			//IJ.log("submit_task(" + i + ")");
			int actual_op = op_mode;
			if(actual_op == 4)
			{	actual_op = 2;} //just like polling only without sf but with the other images
			if( (i >= 0) && (i<=use_threads) /*&& busy_frames[i] >= 0*/)
			{
				//boolean fresh = true;
				if(i < use_threads)
				{	
					++active_frames;
					//fresh = (fr_p_worker[i] == 0);
					++(fr_p_worker[i]);		
				}
				/*else //job.bin
				{
					//job.bin is only written after worker 0 has been fed
					//fresh = (fr_p_worker[0] == 1); 
				}*/
				
				//send all relevant variables
				int job = 0;
				int rank = 0;
				int worker = 0;
				int rID = (i<use_threads) ? busy_frames[i] : monitor_worker;
				int frame = total_reports[rID].frNr;
				int msk = total_reports[rID].mkNr;
				int stab = total_reports[rID].stability - 2*total_reports[rID].runs;
				if(stab < 0) {	stab = 0;}
				int pp = peak_points + total_reports[rID].runs;
				if(pp > 7)	{	pp = 7;}
				if(i < use_threads)
				{
					job = threadVec.get(i).job;
					rank = threadVec.get(i).rank;
					worker = threadVec.get(i).worker + thread_offset;
				}
				//if(fresh) //add extra zero for easier changes in job.bin
				//{
					pwos[i].println("ThreadNum(0" + worker + ")");
					pwos[i].println("#" + Datatitle);
					pwos[i].println("#" + dateFormat.format(date));
					pwos[i].println("ImpWidth(" + impWidth+ ")");
					pwos[i].println("ImpHeight(" + impHeight + ")");
					pwos[i].println("MaskScale(" + mask_scaling + ")");
					pwos[i].println("ModelSize(0"  + modelsize + ")");
					pwos[i].println("SolidSize(0"  + solidsize + ")");
					pwos[i].println("BondLength(0" + bondlength + ")");
					pwos[i].println("MiniLength(0" + minilength + ")");
					pwos[i].println("FieldofView(0" + field_of_view + ")");
					pwos[i].println("MarkSF(" + (mark_sf?1:0) +")");
					pwos[i].println("PeakPoints(0" + pp + ")");
					pwos[i].println("Resolution(0" + uc_res + ")");
					pwos[i].println("Smoothing(0" + smooth_passes + ")");
					pwos[i].println("UseMirror("+ (total_reports[rID].use_mirror ? 1 : 0) + ")");
					pwos[i].println("UsePosition("+ (use_positionQ ? 1 : 0) + ")");
					pwos[i].println("MinHexContrast(" + min_hex_contrast + ")");
					pwos[i].println("RotateTilt("+ (rotate_tilt) + ")");
					pwos[i].println("Task(" + actual_op + ")");
					if(i==use_threads)
					{	
						pwos[i].println("Locking(0)");
						//pwos[i].println("HasConsole(0)");
						pwos[i].println("HasConsole(1)");
						pwos[i].println("ReportUC(0)");
						pwos[i].println("Reporting(1)");
						pwos[i].println("ReportHexImg(0)");
						pwos[i].println("ReportSF(0)");
					}
					else
					{	
						pwos[i].println("Locking(" + (locking?1:0) +")");
						pwos[i].println("HasConsole(0)");
						if( (uc_sum != null) && (uc_sym != null) )
						{	pwos[i].println("ReportUC(1)");}
						else
						{	pwos[i].println("ReportUC(0)");}
						
						if( i == monitor_worker )
						{	pwos[i].println("Reporting(1)");}
						else
						{	pwos[i].println("Reporting(0)");}
						
						if( (hexImg != null) && (sf_sum != null) ) //report hexed frame and sum of hexframes
						{	pwos[i].println("ReportHexImg(1)");} //fairly expensive DEBUG output 
						else
						{	pwos[i].println("ReportHexImg(0)");}
						
						if(	sfImg != null )
						{	pwos[i].println("ReportSF(1)");}
						else
						{	pwos[i].println("ReportSF(0)");}
					}
				//}
				pwos[i].println("#" + total_reports[rID].label);
				pwos[i].println("Frame(" + frame + ")");
				pwos[i].println("Stability(" + stab + ")");
				pwos[i].println("HexImgScale(" + hexImgScale + ")");
				//pwos[i].println("HexEdgeLen("+ total_reports[rID].hexedgelen + ")");
				pwos[i].println("Tilt("+ total_reports[rID].tilt + ")");
				pwos[i].println("EllipseA("+ total_reports[rID].ellipseA + ")");
				pwos[i].println("EllipseB("+ total_reports[rID].ellipseB + ")");
				pwos[i].println("EllipsePhi("+ total_reports[rID].ellipsePhi + ")");
				pwos[i].println("OffsetX("+ total_reports[rID].offset_x + ")");
				pwos[i].println("OffsetY("+ total_reports[rID].offset_y + ")");
				pwos[i].println("OffsetHexQ("+ total_reports[rID].offset_hexq + ")");
				pwos[i].println("OffsetHexR("+ total_reports[rID].offset_hexr + ")");
				pwos[i].println("Min_Light_Dirt(" + min_molecule_coverage + ")");
				pwos[i].println("Max_Light_Dirt(" + max_molecule_coverage + ")");
				pwos[i].println("Light_Dirt_Value(" + molecule_value + ")");
				pwos[i].println("Retune_Origin(" + (retune_origin?1:0) + ")");
				pwos[i].println("Skip_Grid(" + (skip_grid?1:0) + ")");
				pwos[i].println("Create_Defect(0" + create_defect + ")");
				pwos[i].println("NumSF("+ total_reports[rID].subframes + ")");
				//this is also correct
				//System.out.println("Frame: " + frame + "\tFrame: " + total_reports[frame-1].frNr + "\tSubFrames: " + total_reports[frame-1].subframes);
				/*
				pwos[i].println("Merit("+ total_reports[rID].merit + ")");
				pwos[i].println("Raw_Avg("+ total_reports[rID].raw_avg + ")");
				pwos[i].println("Hex_Avg("+ total_reports[rID].hex_avg + ")");
				pwos[i].println("Hex_Low("+ total_reports[rID].hex_low + ")"); 
				pwos[i].println("Contrast("+ total_reports[rID].contrast + ")");
				pwos[i].println("Position_Q("+ total_reports[rID].positionQ + ")");
				pwos[i].println("Mirror_Q("+ total_reports[rID].mirrorQ + ")");
				*/
				if(target_useable && (uc_target != null) )
				{
					pwos[i].println("updateUCTarget(1)");
					send_uc_target(i);
				}
				
				if(has_mask) //send mask before Frame so getFrame in hexsampler can know about it
				{
					pwos[i].println("updateMask(1)");
					sendMask(msk,i);
				}
				pwos[i].println("updateFrameData(1)"); ///maybe the omitted dummy argument crashed the worker?
				sendFrame(frame,i);
				pwos[i].println("Run(1)");
				pwos[i].flush();
				
				if(	i < use_threads)
				{
					String reply = dis[i].readLine();
					while(reply.startsWith("Message") && reader[i].ready())
					{
						IJ.log("worker" + i + ": " + reply.substring(8));
						reply = dis[i].readLine();
					}
					if(!reply.startsWith("LAUNCH"))
					{	
						IJ.log("worker" + i + ": " + reply);
						return false;
					}
					IJ.showStatus("submitted frame " + frame + " to worker" + i );
				}
			}
			else
			{
				System.out.println("illegal request for task submission");
				throw new IllegalStateException("check: '0<=i<=use_threads' FAILED");	
			}
			return true;
		}
		catch(Exception e)
		{
			try
			{
				///There is nothing to read
				if(i<use_threads)
				{
					while(err_reader[i].ready())
					{
						IJ.log( err_reader[i].readLine() );		
					}
				}
			}
			catch(Exception e2)
			{
				IJ.log("failed to read stderr of worker" + i);
				e2.printStackTrace();
			}
			e.printStackTrace();
			return false;
		}
	}

	private boolean setupPolling()
	{ 
		
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog( getClass().getSimpleName() + " " + op_modes[op_mode] );
		gd.setSmartRecording(true);
		gd.addMessage("empty graphene (mask == 1) is always allowed");
		gd.addNumericField("Minimal target coverage", min_molecule_coverage,3);
		gd.addNumericField("Maximal target coverage", max_molecule_coverage,3);
		gd.addNumericField("target value (in mask)", molecule_value,0);
		gd.addNumericField("dump save limit",max_stack_size,0);
		gd.addCheckbox("retune origin", retune_origin || reset_bad_frames);
		gd.addCheckbox("markSF", mark_sf);
		if(!batch_run)
		{
			gd.showDialog();
			if (gd.wasCanceled())
			{	return false;}
		}
		min_molecule_coverage = gd.getNextNumber();
		max_molecule_coverage = gd.getNextNumber();
		molecule_value = (int)gd.getNextNumber();
		max_stack_size = (int)gd.getNextNumber();
		retune_origin = gd.getNextBoolean();
		mark_sf =  gd.getNextBoolean();
		return true;
	}


	private boolean setupDialog()
	{ 
		
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog( getClass().getSimpleName() + " for " + user );
		gd.setSmartRecording(true);
		int[] idArray = WindowManager.getIDList(); // list of all opened images (IDs)
		int idlen = 0;
		if(idArray != null)
		{	idlen = idArray.length;}
		String[] titleArray = new String[idlen + 1]; // titles of opened images and "<open>"
		String[] titleArray2 = new String[idlen+3];  // titles + 
		String[] titleArray3 = new String[idlen+2];  // <open> + <none> + titles
		for (int i = 0; i < idlen; ++i)
		{	
			titleArray3[i+2] = titleArray2[i+1] = titleArray[i] = WindowManager.getImage(idArray[i]).getTitle();		
		}
		titleArray[idlen] = "<open>";
		titleArray2[0] = "<new>";
		titleArray2[idlen+1] = "<open>";
		titleArray2[idlen+2] = "<none>";
		titleArray3[0] = "<none>";
		titleArray3[1] = "<open>";
		if(Datatitle==null)	{	Datatitle = titleArray[0]; }
		
		{
			int fid = 0;
			while(fid < filter_choices.length && !filter_choices[fid].equals(filter_name))
			{	++fid;	}
			if(fid == filter_choices.length) //filter_name is not yet included in filter_choices
			{
				String ofc[] = filter_choices;
				filter_choices = new String[ofc.length + 1];
				for(int f = 0; f < ofc.length; ++f)
				{	filter_choices[f] = ofc[f];}
				filter_choices[ofc.length] = filter_name;
			}
		}
		gd.addMessage(path);
		gd.addChoice("Data input", titleArray, Datatitle);
		titleArray[idlen] = "<new>";
		gd.addChoice("Data Mask", titleArray3, Masktitle);
		gd.addChoice("target unit cell", titleArray2, uc_target_title);
		gd.addNumericField("hexImgScale", hexImgScale, 6);
		gd.addStringField("naming prefix for images", prefix, 30);
		gd.addNumericField("use frames (or entries) (<1 All)", init_use_frames, 0);
		gd.addNumericField("modelsize", modelsize, 0);
		gd.addNumericField("solid (exclusive) size", solidsize, 0);
		gd.addNumericField("Bond Length", bondlength, 0);
		gd.addNumericField("oversampled Bond Length", minilength, 0);
		gd.addNumericField("resolution", uc_res, 0);
		gd.addNumericField("smoothing", smooth_passes, 0);
		gd.addNumericField("stability", stability, 0);
		gd.addNumericField("Num peaks/points", peak_points, 0);
		gd.addNumericField("field_of_view (nm)", field_of_view, 6);
		gd.addNumericField("minimal required contrast for hexagonal mode", min_hex_contrast, 6);
		gd.addChoice("input file",filenames,inputname);
		gd.addChoice("frame filter defs.",filter_choices,filter_name);
		gd.addChoice("task", op_modes, op_modes[op_mode]);
		//gd.addChoice("ruling merit", merits, merits[ruling_merit]);
		//gd.addCheckbox("rearange frames", reorder_stack);
		gd.addCheckbox("skip grid search", skip_grid);
		gd.addCheckbox("reversed start",backwards);
		gd.addCheckbox("sisyphos mode", forever);
		//gd.addCheckbox("use positionQ", use_positionQ);
		//gd.addCheckbox("also use mirrorQ", use_mirrorQ);
		gd.addCheckbox("binary mode", binary_mode);
		gd.addCheckbox("reset states", reset_states);
		
		gd.addNumericField("create defect", create_defect,0);
		gd.addNumericField("rotate unitcells (60° steps)", rotate_tilt,0);
		gd.addNumericField("default hex edge len (<= 0 auto from fov)", def_hel,6);
		gd.addNumericField("default tilt", tilt, 6);
		gd.addNumericField("default excentricity (a/b)", def_excent, 6);
		gd.addNumericField("default phi of ellipse", def_phi,6);
		gd.addNumericField("Number of workers (<1: all)", initthreadnum, 0);
		
		gd.addMessage("0 .. none, 1 .. Warnings & hexed frames, 2 .. also write Job.bin");
		gd.addNumericField("debug level", debug_level, 0);
		gd.addNumericField("monitor worker(-1 none)", monitor_worker, 0);
		
						

		if(!batch_run)
		{
			gd.showDialog();
			if (gd.wasCanceled())
			{	return false;}
		}
		Datatitle = gd.getNextChoice();
		Masktitle = gd.getNextChoice();
		uc_target_title = gd.getNextChoice();
		hexImgScale = gd.getNextNumber();
		prefix = gd.getNextString();
		init_use_frames = (int) gd.getNextNumber();
		modelsize = (int) gd.getNextNumber();
		solidsize = (int) gd.getNextNumber(); 
		bondlength = (int) gd.getNextNumber();
		minilength = (int) gd.getNextNumber();
		uc_res = (int) gd.getNextNumber();
		smooth_passes = (int) gd.getNextNumber();
		stability = (int) gd.getNextNumber(); 
		peak_points = (int) gd.getNextNumber(); 
		field_of_view = gd.getNextNumber();
		min_hex_contrast = gd.getNextNumber();
		inputname = gd.getNextChoice();
		filter_name = gd.getNextChoice();
		String mode = gd.getNextChoice();
		//String ruling_merit_name = gd.getNextChoice();
		//reorder_stack = gd.getNextBoolean();
		
		skip_grid = gd.getNextBoolean();
		backwards = gd.getNextBoolean();
		forever = gd.getNextBoolean();
		//use_positionQ = gd.getNextBoolean();
		//use_mirrorQ = gd.getNextBoolean();
		binary_mode = gd.getNextBoolean();
		reset_states = gd.getNextBoolean();
		
		create_defect = (int)gd.getNextNumber();
		rotate_tilt = (int)gd.getNextNumber();
		def_hel = gd.getNextNumber();
		tilt = gd.getNextNumber();
		def_excent = gd.getNextNumber();
		def_phi = gd.getNextNumber();
		initthreadnum = (int) gd.getNextNumber();
		
		debug_level = (int) gd.getNextNumber(); 
		monitor_worker = (int) gd.getNextNumber();
		
		/*ruling_merit = 0;
		while( !ruling_merit_name.equals(merits[ruling_merit]) )
		{ ++ruling_merit;}
		IJ.log("merit: " + merits[ruling_merit]);*/
		
		op_mode = 0;
		while( !mode.equals(op_modes[op_mode]) )
		{ ++op_mode;}
		IJ.log("task: " + op_modes[op_mode]);
		if( true )
		{
			return setupPolling();
		}
		return true;
	}
	
	private int validateInput()
	{
		
		try
		{
			int res = 0;
			if(Datatitle.equals("<open>")) //duplicate code with fetch_Idata()
			{
				String newdata = IJ.getFilePath("open Idata");
				if(newdata == null)
				{	return 1;}
				if(newdata.endsWith(".zip"))
				{
					Idata = new ImagePlus(newdata);
				}
				else
				{
					Idata = IJ.openVirtual(newdata);
				}
				
				
				Idata.show();
				Thread.sleep(10);
				Datatitle = Idata.getTitle();
				FileInfo fi = Idata.getOriginalFileInfo();
				sourcepath = fi.directory;
				IJ.log("new sourcepath: " + sourcepath);				
			}
			else
			{	
				Idata = WindowManager.getImage(Datatitle);
				FileInfo fi = Idata.getOriginalFileInfo();
				if(fi != null)
				{
					sourcepath = fi.directory;
					IJ.log("new sourcepath: " + sourcepath);
				}	
			}
			if(Idata.getType() != ImagePlus.GRAY16)
			{
				IJ.log(Datatitle + " must be a GRAY16 image!");
				res = 2;
			}
			
			
			if ( (!Masktitle.equals("<none>")) &&
				 (!Masktitle.equals("<open>"))	)
			{	
				Mask = WindowManager.getImage(Masktitle);
				if(Mask != null)
				{
					has_mask = true;
				}
				else
				{
					IJ.log("Could not find Mask " + Masktitle);
					Masktitle = "<open>";
				}
			}
			
			if(Masktitle.equals("<open>")) //duplicate code with fetch_Mask()
			{
				String newdata = IJ.getFilePath("open Mask");
				if(newdata == null)
				{	return 1;}
				if(newdata.endsWith(".zip"))
				{
					Mask = new ImagePlus(newdata);
				}
				else
				{
					Mask = IJ.openVirtual(newdata);
				}
				Mask.show();
				Thread.sleep(10);
				Masktitle = Mask.getTitle();
				has_mask = true;			
			}
			
			has_mask = (Mask!=null);
			
			if(has_mask && (Mask.getType() != ImagePlus.GRAY8) )
			{
				IJ.log(Masktitle + " must be a GRAY8 image");
				res = 2;
			}
			
			
			if( !uc_target_title.startsWith("<") &&
				!uc_target_title.endsWith(">")	)
			{	//a regular Image name
				uc_target = WindowManager.getImage(uc_target_title);
				if(uc_target != null) //it does exist
				{
					target_useable = true;
					update_target = false;
				}
				else
				{
					uc_target_title = "<open>"; //let the user pick a file
				}
			} //no else here since uc_target_title may be changed
			if( uc_target_title.equals("<none>") )
			{
				target_useable = false;
				update_target = false;
				uc_target = null;
			}
			else if ( uc_target_title.equals("<new>") )
			{
				uc_target = NewImage.createShortImage( prefix + "uc_target.tif" , 2*minilength , 2*minilength , 1, NewImage.FILL_BLACK);
				uc_target_title = uc_target.getTitle();
				uc_target.show();
				Thread.sleep(10);
				target_useable = false;
				update_target = true;
			}
			else if(uc_target_title.equals("<open>")) 
			{
				String newtitle = IJ.getFilePath("open unit cell target");
				if(newtitle == null)
				{	return 1;}
				uc_target = new ImagePlus(newtitle);
				uc_target.show();
				Thread.sleep(10);
				uc_target_title = uc_target.getTitle();
				update_target = false;
				target_useable = true;		
			}
			if(uc_target != null)
			{	
				uc_targetSt = uc_target.getStack();
				if(uc_target.getType() != ImagePlus.GRAY16)
				{
					IJ.log(uc_target_title+ " must be a GRAY16 image!");
					res = 2;
				}
				
			}
			
			
			if(inputname.equals("<open>"))
			{
				inputname = IJ.getFilePath("open frame stats file");
				if(inputname == null)
				{	return 1;}
				frame_file = new File(inputname);	
			}
			else if( !( inputname.startsWith("<") && inputname.endsWith(">") ) ) //a regular filename
			{	frame_file = new File(path + inputname );}
			if(frame_file != null)
			{
				boolean ffok = (frame_file.isFile() && frame_file.setReadable(true));
				if(!ffok)
				{
					IJ.log("Error: cannot read " + frame_file.getPath() );
					return 1;
				}
			}
			
			if(filter_name.equals("<open>"))
			{
				filter_name = IJ.getFilePath("open frame filter file");
				if(filter_name == null)
				{	return 1;}
				filter_file = new File(filter_name);	
			}
			else if( !( filter_name.startsWith("<") && filter_name.endsWith(">") ) ) //a regular filename
			{	filter_file = new File(path + filter_name );}
			if(filter_file != null)
			{
				boolean filter_ok = (filter_file.isFile() && filter_file.setReadable(true));
				if(!filter_ok)
				{
					IJ.log("Error: cannot read " + filter_file.getPath() );
					return 1;
				}
			}
			
			if( ( (op_mode == 2) || (op_mode == 3) ) && ( (frame_file == null) || frame_file.getName().equals("frame_init.txt") ) )
			{	
				IJ.log("Error: polling subframes and debug run require any other input file than \"frame_init.txt\"");
				return 1;
			}
			IdataSt = Idata.getStack();
			impWidth = Idata.getWidth();
			impHeight = Idata.getHeight();
			impDepth = IdataSt.getSize();
			if(has_mask)
			{
				MaskSt = Mask.getStack();
				mskDepth = MaskSt.getSize();
				mask_scaling = 1;
				if( (impWidth != Mask.getWidth()) ||
				    (impHeight != Mask.getHeight()) )
				{
					mask_scaling = (int)(impWidth / Mask.getWidth() );
					if( ( impWidth != mask_scaling * Mask.getWidth() ) ||
					    (impHeight != mask_scaling * Mask.getHeight() ) )
					{
						IJ.log("Error: Dimensions of Mask and Data dont match (with scaling " + mask_scaling + ")!");
						return 1;
					}
				}
				 
			}
			
			hexImgScale = Math.abs(hexImgScale);
			if(hexImgScale <= 0.01)
			{	IJ.log("hexImgScale is ridicously small " + hexImgScale);}
			use_frames = init_use_frames;
			if( (use_frames < 1) /*|| (use_frames > impDepth)*/ ) //actually a frame_stats file could contain the same frame with different parameters
			{	use_frames = impDepth;}
			IJ.log("frames/tasks: " + use_frames);
			IJ.log("bondlength: " + bondlength);
			if( (bondlength % 2 != 0) || (bondlength < 1) )
			{
				res = 1; //serious problem
				IJ.log("Error: bondlength must be a positive even integer");
			}
			IJ.log("modelsize: " + modelsize);
			if( (modelsize < 1) || (modelsize % bondlength != 0) )
			{
				res = 1; //serious problem
				IJ.log("Error: modelsize must be a multiple of bondlength");
			}
			IJ.log("solidsize: " + solidsize);
			if( (solidsize < 1) || (solidsize % bondlength != 0) || (solidsize > modelsize) )
			{
				res = 1; //serious problem
				IJ.log("Error: solidsize must be a positive multiple of bondlength and not bigger than modelsize");
			}
			IJ.log("oversampling bondlength: " + minilength);
			///This should actually be possible, let us see
			/*if( (op_mode >= 2) && ( solidsize < modelsize ) && mark_sf )
			{	
				IJ.log("Error: polling AND marking overlapping subframes is not supported");
				return 1;
			}*/
			
			if( (minilength % 2 != 0) || (minilength < 1) )
			{
				res = 1; //serious problem
				IJ.log("Error: oversampling bondlength must be a positive even integer");
			}
			if( minilength < bondlength)
			{
				IJ.log("Warning: sampling bondlength should not be less than actual bondlength");
			}
			if(!uc_target_title.equals("<new>") && !uc_target_title.equals("<none>"))
			{
				if( (uc_target.getWidth() != 2*minilength) ||
					(uc_target.getWidth() != 2*minilength) )
				{
					res = 1;
					IJ.log("Error: " + uc_target_title + " does not match the requested oversampling bondlength");
				}
			}
			if(min_molecule_coverage > max_molecule_coverage)
			{
				res = 1;
				IJ.log("Min_Molecule_Coverage cannot be less than Max_Molecule_Coverage");
			}
			
			if(max_molecule_coverage < 0.0)
			{
				res = 1;
				IJ.log("Max_Molecule_Coverage cannot be negative");
			}
			
			if(field_of_view < 0.0)
			{
				res = 1;
				IJ.log("field_of_view must be positive");
			}
			if(def_hel <= 0.0)
			{	hexedgelen = (impWidth * a0 / field_of_view) / (Math.sqrt(3.0)*bondlength);}
			else
			{	hexedgelen = def_hel;}
			if(def_excent < 1.0)
			{
				res = 1;	
				IJ.log("default excentricity must be >= 1.0");
			}
			def_phi = def_phi % (Math.PI);
			if(def_phi < -Math.PI/2)
			{	def_phi += Math.PI;}
			if(def_phi > Math.PI/2)
			{	def_phi -= Math.PI;}
			
			tilt = tilt % (Math.PI/3);
			if(tilt < -Math.PI/6)
			{	tilt += Math.PI/3;}// stay above -30°
			if(tilt > Math.PI/6)// stay below +30°
			{	tilt -= Math.PI/3;}
			// -30° < rotation <= +30°
			if(peak_points < 1)
			{	peak_points = 1;}
			if(stability < 0)
			{	
				IJ.log("WARNING: minimal supported stability is 0");
				stability = 0;
			}
			if(stability > 3)
			{	
				IJ.log("WARNING: maximal supported stability is 3");
				stability = 3;
			}
			IJ.log("stability: " + stability);
			if(uc_res < 7)
			{	
				uc_res = 7;
				res = 1;
				IJ.log("uc_res cannot be less than 7");
			}
			IJ.log("sampling resolution for unitcells: " + uc_res);
			if(smooth_passes < 0)
			{	
				smooth_passes = 0;
				res = 1;
				IJ.log("Smoothing cannot be less than zero");
			}
			if(smooth_passes > uc_res/2)
			{	
				IJ.log("Warning smoothing exceeds unitcell");
			}
			IJ.log("Smoothing: " + smooth_passes);
			threadnum = initthreadnum;
			if( (threadnum < 1) || (threadnum > max_threads) )
			{	threadnum = max_threads;}
			if( threadnum > use_frames) //we dont need more workers than frames
			{	threadnum = use_frames;}
			IJ.log("workers: " + threadnum);
			if(debug_level < 0 || debug_level > 2)
			{	debug_level = 0;}
			IJ.log("Num peaks/points: " + peak_points);
			IJ.log("debug level: " + debug_level);	
			/*if(((op_mode != 0) && (op_mode != 3)) && forever)
			{
				forever = false;
				IJ.log("cancelled sisyphos mode");
			}*/
			if(forever)
			{	IJ.log("Sisyphos optimization," + (endless?" will never finish": " as long as all workers are busy") );}
			if( (monitor_worker >= 0) && (monitor_worker < threadnum))
			{	IJ.log("console output for worker" + monitor_worker);}
			else
			{	monitor_worker = -1;}
			IJ.log("detect and fix bad frames: " + reset_bad_frames );
			IJ.log("skip (very expensive) grid/peak search: " + skip_grid );
			IJ.log("retune Origin: " + retune_origin);
			if(reset_bad_frames && (!retune_origin))
			{
				IJ.log("WARNING: resetting frames without retuning the origin may be a BAD IDEA");
			}
			IJ.log("use positionQ: " + use_positionQ);
			/*if(use_mirrorQ && !use_positionQ)
			{
				IJ.log("mirrorQ can only be employed when positionQ is in place");
				use_mirrorQ = false;
			}*/
			IJ.log("create defect: " + create_defect);
			IJ.log("use mirrorQ: " + use_mirrorQ);
			IJ.log("binary_mode: " + binary_mode);
			IJ.log("rotate unitcells by 60° steps: " + rotate_tilt);
			if(max_relaunches > 0)
			{
				IJ.log("max_relaunches (per worker): " + max_relaunches);
			}
			if(backwards)
			{
				IJ.log("reversed launch sequence");
			}
			return res;
		}
		catch(Exception e)
		{	
			e.printStackTrace();
		}
		return 1;
	}
	
	private boolean fetch_Idata()
	{
		if( oldDatatitle!=null )
		{ 
			try
			{
			
				if(  ( ! oldDatatitle.equals("<open>") )
				&& ( WindowManager.getImage(oldDatatitle) == null )  )
				{
					File tmp = new File(sourcepath + oldDatatitle);
					if(tmp.isFile())
					{
						IJ.log("opening: " + tmp.getPath());
						if(tmp.getPath().endsWith(".zip"))
						{
							Idata = new ImagePlus(tmp.getPath());
						}
						else
						{
							Idata = IJ.openVirtual(tmp.getPath());
						}
						Idata.show();
						Thread.sleep(10);
					}
					else
					{
						IJ.log("could not find " + tmp.getPath() );
						tmp = new File(sourcepath + oldDatatitle + ".zip");
						if(tmp.isFile())
						{
							IJ.log("opening: " + tmp.getPath());
							Idata = new ImagePlus(tmp.getPath() );
							Idata.show();
							Thread.sleep(10);
						}
						else
						{IJ.log("could not open " + tmp.getPath());}
					}	
				}
				if(	!oldDatatitle.equals("<open>")) 
				{
					if( WindowManager.getImage(oldDatatitle) != null )
					{	Datatitle = oldDatatitle;}
					else  //duplicate code with validateInput()
					{
						if(batch_run) {	IJ.log("could not open Idata during batchrun\n" + oldDatatitle);}
						String newdata = IJ.getFilePath("open Idata");
						if(newdata == null)
						{	return false;}
						Idata = new ImagePlus(newdata);
						Idata.show();
						Thread.sleep(10);
						Datatitle = Idata.getTitle();
						FileInfo fi = Idata.getOriginalFileInfo();
						sourcepath = fi.directory;
						IJ.log("new sourcepath: " + sourcepath);			
					}
				}
			}
			catch(Exception e)
			{	
				e.printStackTrace();
				return false;
			}
			
		}
		return true;
	}
	
	private boolean fetch_Mask()
	{
		if( oldMasktitle!=null )
		{ 
			try
			{
			
				if(  ( ! oldMasktitle.equals("<open>") ) &&
					 ( ! oldMasktitle.equals("<none>") ) &&
				 ( WindowManager.getImage(oldMasktitle) == null )  )
				{
					File tmp = new File(sourcepath + oldMasktitle);
					if(tmp.isFile())
					{
						IJ.log("opening: " + tmp.getPath());
						if(tmp.getPath().endsWith(".zip"))
						{
							Mask = new ImagePlus(tmp.getPath());
						}
						else
						{
							Mask = IJ.openVirtual(tmp.getPath());
						}
						Mask.show();
						Thread.sleep(10);
						has_mask = true;
					}
					else
					{
						IJ.log("could not find " + tmp.getPath() );
						tmp = new File(sourcepath + oldMasktitle + ".zip");
						if(tmp.isFile())
						{
							IJ.log("opening: " + tmp.getPath());
							Mask = new ImagePlus(tmp.getPath() );
							Mask.show();
							Thread.sleep(10);
							has_mask = true;
						}
						else
						{IJ.log("could not open " + tmp.getPath());}
					}	
				}
				if(	( ! oldMasktitle.equals("<open>") ) &&
					( ! oldMasktitle.equals("<none>") )
				  ) 
				{
					if( WindowManager.getImage(oldMasktitle) != null )
					{	Masktitle = oldMasktitle;}
					else  //duplicate code with validateInput()
					{
						if(batch_run) {	IJ.log("could not open Mask during batchrun\n" + oldMasktitle);}
						String newdata = IJ.getFilePath("open Mask");
						if(newdata == null)
						{	return false;}
						Mask = new ImagePlus(newdata);
						Mask.show();
						Thread.sleep(10);
						Masktitle = Mask.getTitle();
						has_mask = true;	
					}
				}
			}
			catch(Exception e)
			{	
				e.printStackTrace();
				return false;
			}
		}
		return true;
	}
	
	private boolean fetch_uc_target()
	{
		if( old_uc_target_title!=null )
		{ 
			try
			{
				
				if(  (  !old_uc_target_title.equals("<open>") && 
						!old_uc_target_title.equals("<new>") && 
						!old_uc_target_title.equals("<none>")    )
				&& ( WindowManager.getImage(old_uc_target_title) == null )  )
				{
					File tmp = new File(path + old_uc_target_title);
					if(tmp.isFile())
					{
						uc_target = new ImagePlus(tmp.getPath());
						uc_target.show();
						Thread.sleep(10);
						update_target = false;
						target_useable = true;
					}
					else
					{
						//IJ.log("could not open " + tmp.getPath());
						tmp = new File(path + old_uc_target_title);
						if(tmp.isFile())
						{
							uc_target = new ImagePlus(tmp.getPath());
							uc_target.show();
							Thread.sleep(10);
							update_target = false;
							target_useable = true;
						}
						else
						{IJ.log("could not open " + tmp.getPath());}
					}	
				}
				if(	!old_uc_target_title.equals("<open>") && 
					!old_uc_target_title.equals("<new>") && 
					!old_uc_target_title.equals("<none>")	) 
				{
					if( (uc_target = WindowManager.getImage(old_uc_target_title)) == null )
					{
						if(batch_run) {	IJ.log("could not open target unit cell during batchrun\n" + old_uc_target_title);}
						String newdata = IJ.getFilePath("open target unit cell");
						if(newdata == null)
						{	return false;}
						uc_target = new ImagePlus(newdata);
						uc_target.show();
						Thread.sleep(10);	
					}
					uc_target_title = uc_target.getTitle();
					update_target = false;
					target_useable = true;	
				}
			}
			catch(Exception e)
			{	
				e.printStackTrace();
				return false;
			}
		}
		return true;
	}
	
	
	
	private boolean readBatchFile()
	{
		IJ.open(batchfile.getPath());
		IJ.selectWindow(batchfile.getName());
		NonBlockingGenericDialog bdg = new NonBlockingGenericDialog("Hex_Magic batchmode");
		bdg.addMessage("You may now edit batchrun.txt and save your changes");
		bdg.addMessage("Click ok to re-read batchrun.txt from disk and proceed");
		bdg.addMessage("The batchrun wont continue once batchrun.txt is closed");
		bdg.addMessage("You may halt and modify a running session with Esc");
		
		bdg.showDialog();
		if(bdg.wasCanceled())
		{	return false;}
		
		try
		{
			Scanner fs = new Scanner(new FileReader(batchfile));
			int repeat = 1;
			boolean skipping = false;
			while( fs.hasNextLine() && !skipping)
			{
				String line = fs.nextLine();
				line = line.trim();
				if(!line.equals("") && !line.startsWith("#"))
				{
					boolean is_a_path = true;
					if(line.startsWith("&"))
					{
						Scanner ls = new Scanner(line);
						String command = ls.next();
						if(command.equals("&REPEAT"))
						{
							if(ls.hasNextInt())
							{	repeat = ls.nextInt();}
							else
							{	IJ.log("Syntax Error in " + batchfile.getPath() + " Integer expected after &REPEAT");}
							is_a_path = false;	
						}
						else if(command.equals("&HERE"))
						{	line = path;}
						else if(command.equals("&EXIT"))
						{
							skipping = true;
							is_a_path = false;
						}
						else
						{	IJ.log("Syntax Error in " + batchfile.getPath() + " unknown command: " + command);}
					}
					if(is_a_path)
					{
						for(int i=0;i<repeat;++i)
						{	batchpaths.add(line);}
						repeat = 1;
						batch_run = true;
					}
				}
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			IJ.showMessage("Error reading " + batchfile.getPath());
			return false;
		}
	
		boolean paths_ok = true;
		for(int batch_id = 0; batch_id < batchpaths.size(); ++batch_id)
		{
			File peek = new File(batchpaths.get(batch_id) + "options.txt");
			if(!peek.isFile())
			{
				IJ.log("missing file: " + peek.getPath());
				paths_ok = false;
			}
		}
		if(!paths_ok)
		{	return false;}
		statusdisplay = new TextWindow("Hex_Sampler status","will run in batch mode as long as batchrun.txt is open",600,300);
		return true;
	}
	
	private boolean loadOptions()
	{
		try
		{
			oldMasktitle = "<none>"; //compatability with pre mask options.txt
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
							{	oldDatatitle = ls.nextLine().trim();}
							else if(token.equals("#Input_Mask"))
							{	oldMasktitle = ls.nextLine().trim();}
							else if(token.equals("#uc_target_title"))
							{	old_uc_target_title = ls.nextLine().trim();}
							else if(token.equals("#filter_name"))
							{	filter_name = ls.nextLine().trim();}
							else if(token.equals("#hexImgScale"))
							{	hexImgScale = ls.nextDouble();}
							else if(token.equals("#bondlength"))
							{	bondlength = ls.nextInt();}
							else if(token.equals("#minilength"))
							{	minilength = ls.nextInt();}
							else if(token.equals("#modelsize"))
							{	modelsize = ls.nextInt();}
							else if(token.equals("#solidsize"))
							{	solidsize = ls.nextInt();}
							else if(token.equals("#peak_points"))
							{	peak_points = ls.nextInt();}
							else if(token.equals("#field_of_view"))
							{	field_of_view = ls.nextDouble();}
							else if(token.equals("#min_hex_contrast"))
							{	min_hex_contrast = ls.nextDouble();}
							else if(token.equals("#min_molecule_coverage"))
							{	min_molecule_coverage = ls.nextDouble();}
							else if(token.equals("#max_molecule_coverage"))
							{	max_molecule_coverage = ls.nextDouble();}
							else if(token.equals("#molecule_value"))
							{	molecule_value = ls.nextInt();}
							else if(token.equals("#retune_origin"))
							{	retune_origin = ls.nextBoolean();}
							else if(token.equals("#create_defect"))
							{	
								if(ls.hasNextBoolean())
								{
									create_defect = ls.nextBoolean()?1:0;
								}
								else
								{
									create_defect = ls.nextInt();
								}
							}
							else if(token.equals("#endless"))
							{	endless = ls.nextBoolean();}
							else if(token.equals("#reset_bad_frames"))
							{	reset_bad_frames = ls.nextBoolean();}
							else if(token.equals("#skip_grid"))
							{	skip_grid = ls.nextBoolean();}
							else if(token.equals("#uc_res"))
							{	uc_res = ls.nextInt();}
							else if(token.equals("#smooth_passes"))
							{	smooth_passes = ls.nextInt();}
							else if(token.equals("#stability"))
							{	stability = ls.nextInt();}
							else if(token.equals("#task"))
							{	String mode = ls.nextLine().trim();
								op_mode = 0;
								while( ( op_mode < op_modes.length ) && !mode.equals(op_modes[op_mode]) )
								{ ++op_mode;}	
								if( op_mode >= op_modes.length )
								{	
									op_mode = 0;
									IJ.log("Error: unknown task: " + mode + "  defaulting to " + op_modes[op_mode]);
								}
							}
							/*else if(token.equals("#ruling_merit"))
							{	String ruling_merit_name = ls.nextLine().trim();
								ruling_merit = 0;
								while( ( ruling_merit < merits.length ) && !ruling_merit_name.equals(merits[ruling_merit]) )
								{ ++ruling_merit;}	
								if( (ruling_merit == merits.length-1) && !ruling_merit_name.equals(merits[ruling_merit]) )
								{	
									ruling_merit = 0;
									IJ.log("Error: unknown ruling_merit: " + ruling_merit_name + "  defaulting to " + merits[ruling_merit]);
								}
							}*/
							else if(token.equals("#hexedgelen"))
							{	def_hel =  ls.nextDouble();}
							else if(token.equals("#tilt"))
							{	tilt = ls.nextDouble();}
							else if(token.equals("#def_excent"))
							{	def_excent = ls.nextDouble();}
							else if(token.equals("#def_phi"))
							{	def_phi = ls.nextDouble();}
							else if(token.equals("#threadnum"))
							{	initthreadnum = ls.nextInt();}
							else if(token.equals("#thread_offset"))
							{	 thread_offset = ls.nextInt();}
							else if(token.equals("#use_frames"))
							{	init_use_frames = ls.nextInt();}
							else if(token.equals("#keep_Idata"))
							{	keep_Idata = ls.nextBoolean();}
							else if(token.equals("#use_mirrorQ"))
							{	use_mirrorQ = ls.nextBoolean();}
							else if(token.equals("#markSF"))
							{	mark_sf = ls.nextBoolean();}
							else if(token.equals("#use_positionQ"))
							{	use_positionQ = ls.nextBoolean();}
							else if(token.equals("#binary_mode"))
							{	binary_mode = ls.nextBoolean();}
							else if(token.equals("#reversed_start"))
							{	backwards = ls.nextBoolean();}
							else if(token.equals("#locking"))
							{	locking = ls.nextBoolean();}
							else if(token.equals("#rotate_tilt"))
							{	
								//FIXME keep these checks as long as there are boolean versions around
								if(ls.hasNextInt()) 
								{	rotate_tilt = ls.nextInt();}
								else if(ls.hasNextBoolean())
								{	rotate_tilt = ls.nextBoolean()?1:0; }
							}
							else if(token.equals("#check_free_nodes"))
							{	check_free_nodes = ls.nextBoolean();}
							else if(token.equals("#debug_level"))
							{	debug_level = ls.nextInt();}
							else if(token.equals("#max_stack_size"))
							{	max_stack_size = ls.nextInt();}
							else
							{
								System.out.println("WARNING unknown token in options.txt: " + token + "\t" + ls.next());
							}
						}
						ls.close();
					}
				}
				fs.close();
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

	void saveOptions(boolean tipps)
	{
		try
		{
			java.io.FileOutputStream opts = new java.io.FileOutputStream(options, false);
			java.io.PrintWriter optp = new java.io.PrintWriter(opts);
			optp.println("options are read from lines that start with #");
			optp.println("Please make sure your filenames do not contain any blank characters");
			if(tipps) optp.println("If images are not yet open they will be read from disk");
			if(tipps) optp.println("Input_Data will only be loaded from sourcepath");
			if(tipps) optp.println("Models and weights from source if they are not found in this directory");
			optp.println("#appname\t" + getClass().getSimpleName());
			optp.println("#sourcepath\t" + sourcepath);
			optp.println("#Input_Data\t" + Idata.getTitle());
			optp.println("#Input_Mask\t" + (has_mask?Mask.getTitle():"<none>"));
			optp.println("#uc_target_title\t" + ( (uc_target!=null) ? uc_target.getTitle() : "<none>" ) );
			optp.println("#filter_name\t" + filter_name);
			optp.println("#hexImgScale\t" + hexImgScale);
			if(tipps) optp.println("stability levels 0 .. tight&risky, 1 .. default, 2 .. wide, 3 .. extensive");
			optp.println("#stability\t" + stability);
			if(tipps) optp.println("the peak is assumed to lie within the bounding box of the N best values");
			optp.println("#peak_points\t" + peak_points);
			if(tipps) optp.println("task to perform on frames with fitting merit");
			optp.println("#task\t" + op_modes[op_mode]);
			//if(tipps) optp.println("decisive merit for optimization");
			//optp.println("#ruling_merit\t" + merits[ruling_merit]);
			if(tipps) optp.println("The bondlength in hexagonal pixels");
			optp.println("#bondlength\t" + bondlength);
			if(tipps) optp.println("The minilength is the resolution used for locating the origin");
			optp.println("#minilength\t" + minilength);
			if(tipps) optp.println("The diameter of hexagonal subframes/models");
			optp.println("#modelsize\t" + modelsize);
			if(tipps) optp.println("The diameter of the solid (exclusive) area of the hexagonal subframes/models");
			optp.println("#solidsize\t" + solidsize);
			if(tipps) optp.println("scale of lateral image dimension in nm");
			optp.println("#field_of_view\t" + field_of_view);
			if(tipps) optp.println("passes for hexagonal next neighbor averaging");
			optp.println("#smooth_passes\t" + smooth_passes);
			if(tipps) optp.println("minimal required contrast to employ hexagonal mode");
			optp.println("min_hex_contrast\t" + min_hex_contrast);
			if(tipps) optp.println("default resolution for sampled unitcells");
			optp.println("#uc_res\t" + uc_res);
			if(tipps) optp.println("hexedgelen <= 0.0 is calculated from FFT and hexedgelen");
			optp.println("#hexedgelen\t" + def_hel);
			if(tipps) optp.println("tilt is between horizontal and armchair direction");
			optp.println("#tilt\t" + tilt);
			if(tipps) optp.println("def_excent is the default excentricity if no frame stats are available");
			optp.println("#def_excent\t" + def_excent);
			if(tipps) optp.println("def_phi is the default angle of the major elliptical axis");
			optp.println("#def_phi\t" + def_phi);
			if(tipps) optp.println("Number of workers to launch for sampling");
			optp.println("#threadnum\t" + initthreadnum);
			if(tipps) optp.println("Number of frames to use from Idata");
			optp.println("#use_frames\t" + init_use_frames);
			if(tipps) optp.println("scan network for available workers, disable it when preparing/launching job files");
			optp.println("#check_free_nodes\t" + check_free_nodes);	
			if(tipps) optp.println("Offset to ID of workers, affects cpu binding");
			optp.println("#thread_offset\t" + thread_offset);
			optp.println("#keep_Idata\t" + keep_Idata);
			if(tipps) optp.println("scale merit with symmmetry / std of uc_pos peaks");
			optp.println("#use_mirrorQ\t" + use_mirrorQ);
			optp.println("#markSF\t" + mark_sf);
			optp.println("#use_positionQ\t" + use_positionQ);
			optp.println("#reversed_start\t" + backwards);
			optp.println("#binary_mode\t" + binary_mode);
			optp.println("#rotate_tilt\t" + rotate_tilt);
			optp.println("#create_defect\t" + create_defect);
			if(tipps) optp.println("if locking is active then workes on the same node quque up for responding to master");
			optp.println("#locking\t" + locking);
			optp.println("#molecule_value\t" + molecule_value);
			optp.println("#min_molecule_coverage\t" + min_molecule_coverage);
			optp.println("#max_molecule_coverage\t" + max_molecule_coverage);
			optp.println("#retune_origin\t" + retune_origin);
			if(tipps) optp.println("endless keeps sisyphus active even if there are less frames than workers");
			optp.println("#endless\t" + endless);
			if(tipps) optp.println("check basic stats after sampling and issue reset if a frame fails");
			optp.println("#reset_bad_frames\t" + reset_bad_frames);
			optp.println("#skip_grid\t" + skip_grid);
			optp.println("debug_levels 0 .. none, 1 .. warnings & hexed frames, 2 .. write job.bin");
			optp.println("#debug_level\t" + debug_level);
			optp.println("#max_stack_size\t" + max_stack_size);
			optp.flush();
			IJ.log("wrote " + options.getPath());
		}
		catch (Exception e)
		{
			System.out.println("Error: could not write to "+ options.getPath());
			e.printStackTrace();
		}
	}

	// returns either a match or the last trial
    private File locateNodes(String[] paths)
    {
		int i = 0;
		File nodes = new File (paths[i] + "nodesSampling.txt");
		while(!(nodes.isFile() && nodes.setReadable(true)) && i < paths.length)
		{
			++i;
			nodes = new File (paths[i] + "nodesSampling.txt");
		}
		return nodes;
	}

	//http://stackoverflow.com/questions/237159/whats-the-best-way-to-check
	//-to-see-if-a-string-represents-an-integer-in-java/237204#237204
	private boolean isInteger(String str) 
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


	int loadNodes()
    {
		String[] places = {path, IJ.getDirectory("imagej")};
		File nodes = locateNodes(places);
		if(nodes.isFile() && nodes.setReadable(true))
		{
			try
			{
				IJ.log(nodes.getPath());
				HashMap<String,Integer> busy_samplers = null;
				boolean has_busy = false;
				
				if(check_free_nodes)
				{
					IJ.log("scanning network for available workers ...");
					File busy_script = new File(IJ.getDirectory("imagej") + "busySamplers.sh");
					if (busy_script.isFile() && busy_script.canExecute())
					{
						has_busy = true;
						busy_samplers = new HashMap<String,Integer>();
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
							if ( !host.equals("<unknown>") && (running != null) && (running.intValue() > 0) && !busy_samplers.containsKey(host) )
							{
								//System.out.println("detected " + running + " workers at " + host);
								busy_samplers.put(host, running);
							}
						}
						ps.destroy(); //also closes all readers
					}
					else
					{
						IJ.log("missing shell script " + busy_script.getPath() );
					}					
				}
				/*else
				{
					IJ.log("Please make sure that these workers are available");
				}*/
				
				threadVec.clear();
				threadadd = 0;
				use_threads = 0;
				
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
					
					if( has_busy && busy_samplers.containsKey(server) && (port >= -1) )
					{
						Integer runInt = busy_samplers.get(server);
						int running = runInt.intValue();
						int more_running = running - multiplicity;
						if( more_running < 0 )
						{	more_running = 0;} 
						busy_samplers.put(server, new Integer(more_running));
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
        //IJ.log("" + multiplicity  + " workers: " + next_command + " performance: " + perf); 
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
		
		for(int i = 0; i < num_threads; ++i)
		{
			threadVec.get(use_threads).job = 0;
			threadVec.get(use_threads).rank = use_threads;
			threadVec.get(use_threads).worker = use_threads;
			++use_threads;
		}
		System.out.println("Workers were pooled");
		return true;
	}
	
	boolean releaseThreads(int num_threads) //ideally called once per job
    {
		if(num_threads > use_threads)
		{
			num_threads = use_threads;
			IJ.log("Error: There are only " + num_threads + " workers to release");
			return false;
		}
		
		if(num_threads < 1)
		{	return false;}
		
		for(int i = 0; (i < num_threads); ++i)
		{
			--use_threads;
			threadVec.get(use_threads).job = -1;
			threadVec.get(use_threads).rank = -1;
			threadVec.get(use_threads).worker = -1;	
		}
		return true;
	}
	
	private boolean saveLog(String path, String logname)
	{
		try
		{
			java.io.FileOutputStream logStream = new java.io.FileOutputStream(path + logname, false);
			java.io.PrintWriter logWriter = new java.io.PrintWriter(logStream);
			//logWriter.println("******* NEW SESSION *******");
			String log = IJ.getLog();
			logWriter.write(log);
			logWriter.close();
			logStream.close();
			return true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return false;
		}
	}
	
		
	void send16bitint(java.io.OutputStream ostr, short val)
    throws java.io.IOException  //low bits first.
    {
        buffer.clear();
        buffer.putShort(0, (short)val);
        ostr.write( bytes, 0, 2);
    }
   
    void send8bitint(java.io.OutputStream ostr, byte val)
    throws java.io.IOException  //low bits first.
    {
        buffer.clear();
        buffer.put(0, (byte)val);
        ostr.write( bytes, 0, 1);
    }
   
    void sendFrame(int fr, int i)
	throws java.io.IOException
	{
		short[] pixels = (short[])IdataSt.getPixels(fr);
		send_short_pixels( i, pixels );
	}

	void sendMask(int mkNr, int i)
	throws java.io.IOException
	{
		byte[] pixels = (byte[])MaskSt.getPixels(mkNr);
		send_byte_pixels( i, pixels );
	}


	void send_uc_target(int i)
	throws java.io.IOException
	{
		uc_targetSt = uc_target.getStack();
		short[] pixels = (short[])uc_targetSt.getPixels(1);
		send_short_pixels( i, pixels );
	}
	
	void send_short_pixels(int i, short[] pixels)
	throws java.io.IOException
	{
		int chsum = 0;
		if(binary_mode)
		{
			pwos[i].println("BeginBinary(" + (2 * pixels.length) + ")");
			pwos[i].flush();	
			for (int j=0; j<pixels.length; ++j)
			{
				int val = (int)(pixels[j] & 0xffff);
				send16bitint(os[i],(short)val);
				chsum += val;
			}
			os[i].flush();
			pwos[i].println("EndBinary(" + chsum + ")");
			pwos[i].flush();
		}
		else //ASCII mode
		{
			pwos[i].println("BeginPixels(" + pixels.length + ")");
			pwos[i].flush();
			for (int j=0; j<pixels.length; ++j)
			{
				int val = (int)(pixels[j] & 0xffff);
				pwos[i].println( Integer.toString(val) );
				chsum += val;
			}
			pwos[i].println("EndPixels(" + chsum + ")");
			pwos[i].flush();
		}
	}
	
	void send_byte_pixels(int i, byte[] pixels)
	throws java.io.IOException
	{
		int chsum = 0;
		if(binary_mode)
		{
			pwos[i].println("BeginBinary(" + (pixels.length) + ")");
			pwos[i].flush();	
			for (int j=0; j<pixels.length; ++j)
			{
				int val = (int)(pixels[j] & 0xffff);
				send8bitint(os[i],(byte)val);
				chsum += val;
			}
			os[i].flush();
			pwos[i].println("EndBinary(" + chsum + ")");
			pwos[i].flush();
		}
		else //ASCII mode
		{
			pwos[i].println("BeginPixels(" + pixels.length + ")");
			pwos[i].flush();
			for (int j=0; j<pixels.length; ++j)
			{
				int val = (int)(pixels[j] & 0xffff);
				pwos[i].println( Integer.toString(val) );
				chsum += val;
			}
			pwos[i].println("EndPixels(" + chsum + ")");
			pwos[i].flush();
		}
	}
	
	
	void update_uc_target()
	{
		if( (uc_target != null) )
		{
			++target_updates;
			int chsum = 0;
			short[] pixels = (short[])uc_targetSt.getPixels(1);
			int sz = uc_miniSt.getSize();
			short[][] ucm_pixels = new short[sz][]; 
			for(int i = 0; i < sz; ++i)
			{
				ucm_pixels[i] = (short[]) uc_miniSt.getPixels(i+1);
			}
			
			for (int j=0; j<pixels.length; ++j)
			{
				double mean = 0.0;
				double weight = 0.0;
				for(int rID=0; rID < use_frames; ++rID)
				{
					double hm = total_reports[rID].hex_merit;
					double dw = Math.pow(0.01*hm,2); //weight with the square of hex_merit
					double el = total_reports[rID].elapsed;
					if( (el > 0.0) && (hm > 0.0) )
					{
						int fr = total_reports[rID].frNr-1;
						weight += dw;
						mean += dw * ucm_pixels[fr][j];
					}	 
				}
				pixels[j] = (weight>0.0)?(short)(mean/weight + 0.5) : 0;
			}
			uc_target.updateAndRepaintWindow();
			try
			{
				Thread.sleep(10);
			}
			catch(Exception e)
			{
				e.printStackTrace();
			}	
			target_needs_update = false;
		}
	}
	
	void dumpsave_sfImg()
	{
			clear_empty_sf(); //There will be empty slices if there has been an earlier dumpsave
			String title = sfImg.getTitle();
			String dumpname = "P" + stack_counter + "_" + title + ".zip";
			
			try
			{
				try
				{
					new FileSaver(sfImg).saveAsZip( path + dumpname );
					System.out.println("saved " + path + dumpname );
				}
				catch(Exception e)
				{
					e.printStackTrace();
					IJ.log("Could not write/create " + path + dumpname);
					return; //sfImg is already harmed
				}
				++stack_counter;
				sfImg.changes = false;
				sfImg.close();
				
				for(int f = 0; f < use_frames; ++f)
				{
					int sl = total_reports[f].frNr;
					String sl_lbl = IdataSt.getSliceLabel(sl);
					sfImgArray[f].changes = false;
					sfImgArray[f].close();
					sfImgArray[f] = NewImage.createShortImage( prefix + sl_lbl, modelsize , modelsize , 1, NewImage.FILL_BLACK);
				}
				sfImg = NewImage.createShortImage( title, modelsize , modelsize , 1, NewImage.FILL_BLACK);
				sfImg.show();
				sfImgSt = sfImg.getStack();
				sfImg.updateAndRepaintWindow();
				Thread.sleep(10);
			}
			catch(Exception e)
			{
				e.printStackTrace();
				IJ.log("Exception while reseting " + title);
			}
			return;
	}


}
