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



public class Hex_Collector implements PlugIn 
{
	
	static java.text.DecimalFormat format = new DecimalFormat("0.000000");
    static DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	static Date date = new Date();
	
	FrameFilter filter = null;
	
	private class Sfbin
	{
		Sfbin() {};
		ImagePlus imp = null;
		ImageStack impSt = null;
		ImagePlus probMap = null;
		ImagePlus model = null;
		ImagePlus asymodel = null;
		boolean first_slice = true;
		double lower_hex_std = -1;
		double upper_hex_std = -1;
		double lower_hex_avg = -1;
		double upper_hex_avg = -1;
		
		boolean accept_frame(int fr)
		{
			final int repID = repIDs[fr-1];
			if(repID == -1 )
			{
				return false;
			}
			HexSampleReport rep = total_reports[ repID ];
			final double hex_std = rep.hex_std; // _scaled
			final double hex_avg = rep.hex_avg;
			boolean hex_std_fit = ( ( (lower_hex_std == -1) || (lower_hex_std <= hex_std ) ) && // _scaled
								   ( (upper_hex_std == -1) || (upper_hex_std > hex_std ) ) ); // _scaled
			boolean hex_avg_fit = ( ( (lower_hex_avg == -1) || (lower_hex_avg <= hex_avg ) ) && // _scaled
								   ( (upper_hex_avg == -1) || (upper_hex_avg > hex_avg ) ) ); // _scaled
			
			
			boolean contrast_ok = (rep.contrast >= min_contrast);
			
			return (hex_std_fit && hex_avg_fit && contrast_ok);					   
		}
		
		void pick_frame(int fr)
		{
			short[] pixels = (short[])IdataSt.getPixels(fr);
			String lbl = IdataSt.getSliceLabel(fr);
			//IJ.log(lbl + " -> " + imp.getTitle());
			if(first_slice)
			{
				impSt.setPixels(pixels,1);
				impSt.setSliceLabel(lbl,1);
				//imp.setStack(impSt);
				first_slice = false;
			}
			else
			{	
				//check if there is already an older version
				int lastslice = impSt.getSize();
				int sl = (check_duplicates?1:(lastslice+1) );
				while(  ( sl <= lastslice ) && ( !impSt.getSliceLabel(sl).equals(lbl) )  )
				{	++sl;}
				if(sl > lastslice) //a new slice
				{	impSt.addSlice(lbl,pixels);}
				else // an update for an earlier slice
				{
					impSt.setPixels(pixels,sl);
					impSt.setSliceLabel(lbl,sl);
				}
				//imp.setStack(impSt);	
			}
		}	
	}
	
	
	ImagePlus tmpImp = null; //Hack for passing references 
	String Datatitle = null;
	ImagePlus Idata = null;
	ImageStack IdataSt = null;
	
	String Othertitle = "uncategorized.tif";
	ImagePlus Odata = null;
	ImageStack OdataSt = null;
	
	String Filteredtitle = "filtered.tif";
	ImagePlus Fdata = null;
	ImageStack FdataSt = null;
	
	//boolean full_octave = false;
	//boolean apply_scaling = false; //no good idea, gives weired Histogramms
	boolean interpolate = true;
	boolean keep_step = true;
	boolean sym_Histo = true;
	
	Sfbin[] sfbins = null;
	String[] probArray = {"Poisson","PoissonF","GammaF","GaussF","Exp","Histo","invGammaF","invBetaF","Rubber","LogNormF"};
	String PDF = "GammaF";
	String frame_name = "<open>";
	String bin_name = "<open>";
	String user = System.getProperty("user.name");
	String prefix = "";
	File options = null;
	File frame_file = null;
	File filter_file = null;
	File bin_file = null;
	HexSampleReport[] total_reports = null;
	String path = IJ.getDirectory("imagej");
	String sourcepath = path;
	String filter_name = "<none>";
	String[] filter_choices = {"<none>","<open>","filter_frames.txt"};
	
	boolean do_dialog = true;
	boolean keep_Idata = false;
	boolean unique_binning = true;
	boolean check_duplicates = false;
	int new_avg_level = 128;
	int new_max_level = 512;
	double min_hex_std;
	double max_hex_std;
	double min_hex_avg;
	double max_hex_avg;
	int bondlength = 4;
	int impDepth = 0;
	int impWidth = 0;
	int impHeight = 0;
	int min_histwidth = 0;
	int[] repIDs = null;
	int modelGraylevels = 128; //reasonable default
	int minfactor = 3;//smallest believed histogram count 
	double modelBrightness = -4; //good default for blank graphene
	double pvalfloor = -3.0; //Roughly 10^-6 minimum prob
	double min_contrast = 0.005;
	double gb_sigma = 5;
	double divide_raw = 1;
	boolean[] hex_mask = null;
	
	int cx = 0;
	int cz = 0;
	int cy = 9;
	
	boolean check_filter(int fr)
	{
		int repID = repIDs[fr-1];
		if(repID == -1 )
		{
			return false;
		}
		
		if(filter == null)
		{	return true;} //trivial case
		HexSampleReport rep = total_reports[ repID ];
		return filter.apply_filter(rep); //only true if all parameters are within the bounds				   
	}
	
	
	private void setHexcenter(int q, int r)
	{
		//in principle negative centers are allowed
		cx = q - (r - (r&1))/ 2;
		cz = r;
		cy = -cx - cz;
	}
	
	
	
	private boolean iswithinHex(int q, int r, int cradius)
	{
		//if(q < 0 || r < 0)
		//{	IJ.error("Hex_Histogramm.iswithinHex: image coordinates cannot be negative");}
		//cube coords of target hex
		int px = q - (r - (r&1)) / 2;
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
	
	
	private void make_hex_mask(int rH)
	{
		setHexcenter(rH, rH);
		hex_mask = new boolean[4 * rH * rH];
		for(int ind = 0; ind < hex_mask.length ; ++ind)
		{	
			int q = ind%(2*rH);
			int r = ind/(2*rH);		
			hex_mask[ind] = iswithinHex(q,r,rH);		
		
		}
	}
	
	
	public void run(String arg)
    {
		IJ.register(getClass());
		Locale.setDefault(Locale.UK);
		date.setTime(System.currentTimeMillis());
        IJ.log(getClass().getSimpleName() + " " + dateFormat.format(date) + "  " + user);
		path = IJ.getDirectory( getClass().getSimpleName() + " for " + user);
		if(path == null) {	return;}
        IJ.log(path);
        options = new File(path + "options.txt");
        if( !loadOptions() )
        {	return;}
        tmpImp = Idata;
        if (!fetch_image(Datatitle))
        {	Datatitle = null;}
		Idata = tmpImp;
		if(Idata != null)
        {	Idata.show();}
		int res = 0;
		do
		{
			if(!setupDialog())	{	return;}
			res = validateInput();
			if(res > 0)	{	IJ.error("bad input! Better luck next time");}
		} while (res > 0);
		
		System.out.println(getClass().getSimpleName() + "\t*********NEW SESSION*********");
		saveOptions(false);
		IJ.log("saved " + options.getPath());
		date.setTime(System.currentTimeMillis());
		IJ.log( "launching: " + dateFormat.format(date));
		IJ.showProgress(0.0);
		IJ.showStatus( getClass().getSimpleName() );
		apply_mean_shift();
		distribute_frames();
		IJ.showProgress(1.0);
		createPtables();
		IJ.showStatus(getClass().getSimpleName() + " done");
		date.setTime(System.currentTimeMillis());
		IJ.log( "finished: " + dateFormat.format(date));
		saveAllBins(true);//ask_user
		System.out.println( getClass().getSimpleName() + "\t*********END SESSION*********");
        
	}

	private boolean createPtables()
	{
		HexAnalyzer ha = new HexAnalyzer();
		if(min_histwidth > 0)
		{	ha.min_histwidth = min_histwidth;}
		int biggest_bin = 0;
		int max_size = 0;
		for(int sfb = 0; sfb < sfbins.length; ++sfb)
		{
			int sz = sfbins[sfb].impSt.getSize();
			if(sz > max_size)
			{
				max_size = sz;
				biggest_bin = sfb;
			}
			IJ.log(sfbins[sfb].imp.getTitle() + "  slices: " + sz);
		}
		
		for(int sfb = 0; sfb < sfbins.length; ++sfb)
		{
			String lbl = sfbins[sfb].imp.getTitle();
			String prob_map_title = prefix + "Prob_Map_" + lbl.substring(prefix.length());
			String model_title = prefix + "Model_" + lbl.substring(prefix.length());
			String asym_model_title = prefix + "AsyModel_" + lbl.substring(prefix.length());
			sfbins[sfb].probMap = ha.get_prob_map(	lbl, PDF, prob_map_title, sym_Histo, 0, bondlength, modelGraylevels,
													minfactor, modelBrightness, pvalfloor, gb_sigma, interpolate, keep_step);
			sfbins[sfb].probMap.show();
			if(sfb == biggest_bin)
			{
				//relies on an earlier call of HexAnalyzer::get_prob_map
				sfbins[sfb].model = ha.make_models(model_title,1); 
				sfbins[sfb].model.show();
				
				sfbins[sfb].asymodel = ha.make_asym_models(asym_model_title,1);
				sfbins[sfb].asymodel.show(); 
			}
		}
		return true;
	}
	
	private boolean saveAllBins(boolean ask_user)
	{
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(getClass().getSimpleName() + " for " + user);
		gd.setSmartRecording(true);
		gd.addMessage(path);
		gd.addMessage("Press Ok to save all non empty stacks and the Prob_Maps now\n" + 
					"They will remain open in any case");
		
		if(ask_user)
		{
			gd.showDialog();
			if (gd.wasCanceled())
			{	return false;}
		}
		if( (!ask_user) || gd.wasOKed())
		{
			System.out.println(path);
			for(int sfb = 0; sfb < sfbins.length; ++sfb)
			{
				int binsize = sfbins[sfb].imp.getStackSize();
				if( sfbins[sfb].imp.getStackSize() > 1 )
				{
					
					saveImg(sfbins[sfb].imp, path);
					saveImg(sfbins[sfb].probMap, path);
					if(sfbins[sfb].model != null)
					{	saveImg(sfbins[sfb].model, path);}
					if(sfbins[sfb].asymodel != null)
					{	saveImg(sfbins[sfb].asymodel, path);}
				}
				else
				{
					sfbins[sfb].imp.close();
					sfbins[sfb].probMap.close();
					if(sfbins[sfb].model != null)
					{	sfbins[sfb].model.close();}
					if(sfbins[sfb].asymodel != null)
					{	sfbins[sfb].asymodel.close();}
				}
			}
		}
		return true;
	}
	
	void saveImg(ImagePlus imgp, String path)
    {
		if(imgp.getStackSize() > 1)
		{	new FileSaver(imgp).saveAsTiffStack( path + imgp.getTitle() );}
		else
		{	new FileSaver(imgp).saveAsTiff( path + imgp.getTitle() );}
		System.out.println("saved " + imgp.getTitle() );
	}
	
	
	
	private void debug_sf_identification()
	{
		int lastID = -1;
		int sf_count = 0;
		int unknown = 0;
		for( int f = 0; f < repIDs.length; ++f)
		{
			int repID = repIDs[f];
			if(repID != -1)
			{
				if( (repID != lastID) || (f == repIDs.length-1) )
				{
					if(f != 0)
					{	IJ.log( "multiplicity: " + sf_count);}
					IJ.log(	total_reports[ repID ].label );
					lastID = repID;
					sf_count = 1;
				}
				else
				{	++sf_count;}
			}
			else
			{	++unknown;}
		}
		IJ.log("unknown subframes: " + unknown);
	}
	
	private void apply_mean_shift()
	{
		if((new_avg_level < 0) && (new_max_level < 0))
		{	return;} // nothing to do;
		if(keep_Idata)
		{
			ImagePlus IdataS = Idata.duplicate();
			Idata = IdataS;
			Idata.show();
			IdataSt = Idata.getStack();
		}
		Idata.setTitle("shifted_" + prefix + Datatitle);
		int radius = impWidth/2;
		long underflow = 0;
		long regular = 0;
		long chopped = 0;
		//long scaled = 0;
		//TODO prefetch hex mask here
		make_hex_mask(radius);
		
		for(int sl = 1; sl <= impDepth; ++sl)
		{
			int repID = repIDs[sl-1];
			if( (repID >= 0) && (repID < total_reports.length) )
			{
				double hex_avg = total_reports[ repID ].hex_avg;
				double hex_std = total_reports[ repID ].hex_std;
				double contr   = total_reports[ repID ].contrast; 
				if( contr >= min_contrast)
				{
					/*
					double scaling = 1.0;
					if(apply_scaling)
					{
						if(hex_std > max_hex_std)
						{
							double inv_scaling = Math.ceil(hex_std/max_hex_std);
							scaling = 1.0/inv_scaling;
							++scaled;
						}
						else if (hex_std < min_hex_std)
						{
							scaling = Math.ceil(min_hex_std/hex_std);
							++scaled;
						}
					}
					*/  
					//total_reports[ repID ].hex_std_scaled = total_reports[ repID ].hex_std; // * scaling 
					int avg_shift = (int)( (new_avg_level > 0) ? (new_avg_level) : hex_avg );
				
				
					//if( (avg_shift != 0) || (scaling != 1.0) ) //actually always true in real world use case
					//{
					//IJ.log("slice: " + sl + " hex_std: " + hex_std + " scaling: " + scaling);
					
					short[] pixels = (short[])IdataSt.getPixels(sl);
					for(int ind = 0; ind < pixels.length; ++ind)
					{
						if(hex_mask[ind]) //relies on earlier call with correct radius
						{
							int q = ind % impWidth;
							int r = ind / impWidth;
							//if(Hex_Magic.iswithinHex(q, r, radius)) //TODO use hex_mask here
							
							pixels[ind] = (short)(((double)pixels[ind] - hex_avg)/divide_raw + avg_shift + 0.5); // * scaling
							
							if(pixels[ind] >= 0)
							{	
								++regular; 	
							}
							else
							{	
								++underflow;
								pixels[ind] = 0;
							}
							if( (new_max_level > 0) && (pixels[ind] > (short)new_max_level))
							{
								++chopped;
								pixels[ind] = (short)new_max_level;
							}
						}
					}
					//}
				}
			}
			
		}
		if(underflow > 0)
		{
			IJ.log("WARNING: zeroed " + underflow + " of " + (regular+underflow) + " pixels in " + Idata.getTitle() );
		}
		if(chopped > 0)
		{
			IJ.log("WARNING: chopped " + chopped + " of " + (regular+underflow) + " pixels in " + Idata.getTitle() );
		}
		/*
		if(scaled > 0)
		{
			IJ.log("applied scaling to " + (scaled) + " of " + impDepth + " slices");
		}*/
		Idata.updateAndRepaintWindow();
	}
	
	
	private boolean setupDialog()
	{ 
		
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(getClass().getSimpleName() + " for " + user);
		gd.setSmartRecording(true);
		int[] idArray = WindowManager.getIDList(); // list of all opened images (IDs)
		int idlen = 0;
		if(idArray != null)
		{	idlen = idArray.length;}
		String[] titleArray = new String[idlen + 1]; // titles of opened images
		for (int i = 0; i < idlen; ++i)
		{	titleArray[i] = WindowManager.getImage(idArray[i]).getTitle();}
		titleArray[idlen] = "<open>";
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
		gd.addStringField("frame stats file", frame_name, 30);
		gd.addChoice("frame filter defs.",filter_choices,filter_name);
		gd.addMessage("global manipulations prior to binning");
		gd.addNumericField("divide raw data", divide_raw,4);
		gd.addNumericField("Shift avg to (<0 no shift)", new_avg_level, 0);
		gd.addNumericField("Chop values at (<0 no limit)", new_max_level, 0);
		gd.addNumericField("Minimal required contrast", min_contrast, 3);
		gd.addNumericField("Fixed range for Prob. Maps (0 auto)", min_histwidth,0);
		gd.addCheckbox("keep original Idata", keep_Idata);
		gd.addStringField("prefix for new images", prefix, 30);
		gd.addStringField("binning definitions (or <open>)", bin_name, 30);
		gd.addCheckbox("validate unique_binning", unique_binning);
		gd.addCheckbox("scan for duplicate subframes(aka updating)", check_duplicates);
		/// Avery bad idea for Gamma Distributions, it may thow be applicable for high dose Gauss distributions
		//gd.addCheckbox("scale data to range (requries full octave)", apply_scaling );
		gd.addStringField("uncategorized frames", Othertitle, 30);
		gd.addNumericField("bondlength", bondlength, 0);
		gd.addNumericField("Model gray levels: ", modelGraylevels, 0);
		gd.addNumericField("Equivalent max data in model1(x<0 x*mean)", modelBrightness, 3);
		gd.addChoice("Noise Model ", probArray ,PDF);
		gd.addCheckbox("symmetrized Histogramm", sym_Histo);
		gd.addCheckbox("interpolate ", interpolate);
		gd.addCheckbox("keep step", keep_step);
		gd.addNumericField("threshold value for histogram", minfactor,0);
		gd.addNumericField("Gaussian Blur sigma (<0 none)", gb_sigma, 3);
		gd.addNumericField("regularization for log probs",pvalfloor, 3);
		if( do_dialog)
		{
			gd.showDialog();
			if (gd.wasCanceled())
			{	return false;}
		}
		Datatitle = gd.getNextChoice();
		frame_name = gd.getNextString();
		filter_name = gd.getNextChoice();
		divide_raw = gd.getNextNumber();
		new_avg_level = (int)gd.getNextNumber();
		new_max_level = (int)gd.getNextNumber();
		min_contrast = gd.getNextNumber();
		min_histwidth = (int)gd.getNextNumber();
		keep_Idata = gd.getNextBoolean();
		prefix = gd.getNextString();
		bin_name = gd.getNextString();
		unique_binning = gd.getNextBoolean();
		check_duplicates = gd.getNextBoolean();
		//apply_scaling = gd.getNextBoolean();
		Othertitle = gd.getNextString();
		gd.addMessage("Parameters for ptables");
		bondlength = (int)gd.getNextNumber();
		modelGraylevels = (int) gd.getNextNumber();
		modelBrightness = gd.getNextNumber();
		PDF = gd.getNextChoice();
		sym_Histo = gd.getNextBoolean();
		interpolate = gd.getNextBoolean();
		minfactor = (int)gd.getNextNumber();
		gb_sigma = gd.getNextNumber();
		pvalfloor = gd.getNextNumber();
		return true;
	}
	
	private int validateInput()
	{
		int res = 0;
		if(Datatitle.equals("<open>"))
		{
			String newdata = IJ.getFilePath("open Idata");
			if(newdata == null)
			{	return 1;}
			Idata = new ImagePlus(newdata);
			Idata.show();
			Datatitle = Idata.getTitle();
			FileInfo fi = Idata.getOriginalFileInfo();
			sourcepath = fi.directory;
			IJ.log("new sourcepath: " + sourcepath);				
		}
		else
		{	Idata = WindowManager.getImage(Datatitle);}
		
		IdataSt = Idata.getStack();
		impWidth = Idata.getWidth();
		impHeight = Idata.getHeight();
		impDepth = IdataSt.getSize();
		if(impWidth != impHeight)
		{
			IJ.log("Error: " + Datatitle + " ought to be square shaped");
			return 1;
		}
		Hex_Magic.setHexcenter(impWidth/2,impHeight/2);
		
		
		if(frame_name.equals("<open>"))
		{
			frame_name = IJ.getFilePath("open frame stats file");
			if(frame_name == null)
			{	return 1;}
			frame_file = new File(frame_name);	
		}
		else if( !( frame_name.startsWith("<") && frame_name.endsWith(">") ) ) //a regular filename
		{	frame_file = new File(sourcepath + frame_name );}
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
			filter_name = filter_file.getName();	
		}
		else if( !( filter_name.startsWith("<") && filter_name.endsWith(">") ) ) //a regular filename
		{	filter_file = new File(sourcepath + filter_name );}
		if(filter_file != null)
		{
			boolean filter_ok = (filter_file.isFile() && filter_file.setReadable(true));
			if(!filter_ok)
			{
				IJ.log("Error: cannot read " + filter_file.getPath() );
				return 1;
			}
			filter = new FrameFilter(filter_file);
		}
		
		int good_sf = loadReports();
		if(good_sf == -1)
		{
			IJ.log("Error while reading " + frame_name);
			return 1;
		}
		
		if(impDepth != good_sf)
		{	IJ.log("Warning: Not all subframes in " + Datatitle + " are descriped in " + frame_name);}
		
		if(bin_name.equals("<open>"))
		{
			bin_name = IJ.getFilePath("open binning description file");
			if(bin_name == null)
			{	return 1;}
			bin_file = new File(bin_name);	
		}
		else if( !( bin_name.startsWith("<") && bin_name.endsWith(">") ) ) //a regular filename
		{	bin_file = new File(sourcepath + bin_name );}
		if(bin_file != null)
		{
			boolean ffok = (bin_file.isFile() && bin_file.setReadable(true));
			if(!ffok)
			{
				IJ.log("Error: cannot read " + bin_file.getPath() );
				return 1;
			}
		}
		
		int good_bins = loadbins();
		if(	good_bins == -1)
		{
			IJ.log("Error while reading " + bin_name);
			return 1;
		}
		if(	good_bins == 0)
		{
			IJ.log("Warning there are no bin_stacks defined in " + bin_name);
		}
		
		if(minfactor < 0)
		{
			IJ.log("Error threshold for interpreting histogram as probability must be >= 0 and not " + minfactor);
			return 1;
		}
		else if (minfactor == 0)
		{
			IJ.log("Warning: zeros in histogram will not be thresholded. Brace for NaNs and Infs!");
		}
		
		
		if(new_avg_level >= 0)
		{	IJ.log("shifting mean to " + new_avg_level);}
		if(new_max_level > 0)
		{	IJ.log("chopping pixel values at " + new_max_level);}
		
		return res;
	}
	
	private int loadReports()
	{
		try
		{
			int numSF = 0;
			repIDs = new int[impDepth]; //every subframe belongs to one sampling report
			for(int i = 0; i < impDepth; ++i)
			{	repIDs[i] = -1;} //mark everyone as not yet identified
			Vector<HexSampleReport> reportVec = new Vector<HexSampleReport>();
			Scanner ffs = new Scanner(new FileReader(frame_file));
			IJ.log("reading " + frame_file.getPath());
			int repID=-1;
			while ( ffs.hasNextLine() )
			{
				String line = ffs.nextLine();
				line = line.trim();
				while (ffs.hasNextLine() && (line.startsWith("#") || line.length() == 0) )
				{
					line = ffs.nextLine();
					line = line.trim();
					if(line.startsWith("#label"))
					{	HexSampleReport.set_key_line(line.substring(1));}
				}
				++repID;
				HexSampleReport report = new HexSampleReport();
				reportVec.add(report);
				report.parseString(line,Idata,null);
				String lbl = report.label;
				int first = report.frNr;
				int last = report.frNr + report.subframes;
				for (int f = first; f < last; ++f)
				{
					String sl = IdataSt.getSliceLabel(f);
					if( sl.contains(lbl))
					{	
						repIDs[f-1] = repID;
						++numSF;	
					}
				}
			}
			++repID;
			total_reports = new HexSampleReport[repID];
			for(int id = 0; id < repID; ++id)
			{	total_reports[id] = (HexSampleReport)reportVec.get(id); }
			//debug_sf_identification();
			return numSF;
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return -1;
		}
	}
	
	private int loadbins()
	{
		try
		{
			Vector<Sfbin> sfbinVec = new Vector<Sfbin>();
			Scanner ffs = new Scanner(new FileReader(bin_file));
			IJ.log("reading " + bin_file.getPath());
			int binID=0;
			
			Sfbin sfbin = null;
			while ( ffs.hasNextLine() )
			{
				String line = ffs.nextLine();
				line = line.trim();
				while (ffs.hasNextLine() && (!line.startsWith("#") || line.length() == 0) )
				{
					//System.out.println("ignore: " + line);
					line = ffs.nextLine();
					line = line.trim();
				}
				//System.out.println("valid: " + line);
				Scanner ls = new Scanner(line);
				if(ls.hasNext())
				{
					String token = ls.next();
					token = token.trim();
					if (token.startsWith("#"))
					{
						if(token.equals("#sub_stack"))
						{
							String tag = ls.nextLine().trim();
							++binID;
							sfbin = new Sfbin();
							sfbinVec.add(sfbin);
							tmpImp = sfbin.imp;
							boolean impfound = fetch_image( prefix + tag);
							sfbin.imp = tmpImp;
							if(impfound)
							{	
								
								if( (sfbin.imp.getWidth() == impWidth) &&
								    (sfbin.imp.getHeight()== impHeight) )
								{	sfbin.first_slice = false;}
								else
								{
									IJ.log("Warning: Dimensions of preexisting " + prefix + tag + " dont match");
									sfbin.imp.setTitle("MISMATCH " + prefix + tag);
									sfbin.imp.show();
									sfbin.imp = null; 
								}
							}
							
							if(sfbin.first_slice)
							{
								sfbin.imp = NewImage.createShortImage( prefix + tag  , impWidth , impHeight , 1, NewImage.FILL_BLACK);
							}
							sfbin.impSt = sfbin.imp.getStack();
							sfbin.imp.show();
						}
						else if (sfbin != null)
						{
							if(	token.equals("#lower_hex_std") )
							{	sfbin.lower_hex_std = ls.nextDouble();}
							else if(	token.equals("#upper_hex_std") )
							{	sfbin.upper_hex_std = ls.nextDouble();}
							else if(	token.equals("#lower_hex_avg") )
							{	sfbin.lower_hex_avg = ls.nextDouble();}
							else if(	token.equals("#upper_hex_avg") )
							{	sfbin.upper_hex_avg = ls.nextDouble();}
							else
							{	IJ.log("Unknown token " + token + " in " + bin_name);}
						}
						else
						{
							IJ.log("Error in " + bin_name + " " + token + " is not preceeded by any #sub_stack");
						}
						
					}
					else
					{
						IJ.log("WTF: " + line);
					}
				}	
			}
			sfbins = new Sfbin[binID];
			boolean first_bin = true;
			for(int id = 0; id < binID; ++id)
			{	
				sfbins[id] = (Sfbin)sfbinVec.get(id);
				if(first_bin)
				{
					min_hex_std = sfbins[id].lower_hex_std;
					max_hex_std = sfbins[id].upper_hex_std;
					min_hex_avg = sfbins[id].lower_hex_avg;
					max_hex_avg = sfbins[id].upper_hex_avg;
					first_bin = false;
				}
				
				
				
				if( sfbins[id].lower_hex_std < min_hex_std )
				{
					min_hex_std = sfbins[id].lower_hex_std;
				}
				if( sfbins[id].upper_hex_std > max_hex_std )
				{
					max_hex_std = sfbins[id].upper_hex_std;
				}
				
				if( sfbins[id].lower_hex_avg < min_hex_avg )
				{
					min_hex_avg = sfbins[id].lower_hex_avg;
				}
				if( sfbins[id].upper_hex_avg > max_hex_avg )
				{
					max_hex_avg = sfbins[id].upper_hex_avg;
				}		
			
			
			
			}
			//full_octave = (max_hex_std == 2* min_hex_std);
			Odata = NewImage.createShortImage( prefix + Othertitle  , impWidth , impHeight , 1, NewImage.FILL_BLACK);
			OdataSt = Odata.getStack();
			Odata.show();
			Fdata = NewImage.createShortImage( prefix + Filteredtitle  , impWidth , impHeight , 1, NewImage.FILL_BLACK);
			FdataSt = Fdata.getStack();
			Fdata.show();
			IJ.log("minimal accepted hex std: " + min_hex_std + " maximal hex std: " + max_hex_std ); //+ " full octave: " + full_octave);
			IJ.log("minimal accepted hex avg: " + min_hex_avg + " maximal hex avg: " + max_hex_avg );
			/*if(apply_scaling)
			{
				apply_scaling &= (full_octave && (new_avg_level >= 0) );
				IJ.log("scaling data to match hex_std: " + apply_scaling);
			}*/
			
			return binID;
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return -1;
		}
	}
	
	private boolean distribute_frames()
	{
		try
		{
			boolean first_other_slice = true;
			boolean first_filtered_slice = true;
			int progStep = impDepth/100;
			if(progStep == 0) { progStep = 1;}
			for(int fr = 1; fr <= impDepth; ++fr)
			{
				boolean taken = false;
				boolean filter_ok = check_filter(fr);
				
				if(filter_ok)
				{
					for(int sfb = 0; sfb < sfbins.length; ++sfb)
					{
						Sfbin sfbin = sfbins[sfb];
						
						boolean take = sfbin.accept_frame(fr);
						if(take)
						{
							if(unique_binning)
							{
								if (taken)
								{	IJ.log("Warning ambiguous binning frame " + fr + ": " + IdataSt.getSliceLabel(fr) );} 
								else
								{	
									sfbin.pick_frame(fr);
									taken = true;
								}
							}
							else
							{	
								sfbin.pick_frame(fr);
								taken = true;
							}		
						}
					}
				}
				if(!taken)
				{
					short[] pixels = (short[])IdataSt.getPixels(fr);
					String lbl = IdataSt.getSliceLabel(fr);
					if(filter_ok)
					{
						if(first_other_slice)
						{
							OdataSt.setPixels(pixels,1);
							OdataSt.setSliceLabel(lbl,1);
							first_other_slice = false;
						}
						else
						{	
							OdataSt.addSlice(lbl,pixels);
						}
					}
					else
					{
						if(first_filtered_slice)
						{
							FdataSt.setPixels(pixels,1);
							FdataSt.setSliceLabel(lbl,1);
							first_filtered_slice = false;
						}
						else
						{	
							FdataSt.addSlice(lbl,pixels);	
						}
					}
				}
				if(fr%progStep == 0)
				{
					IJ.showProgress(fr,impDepth);
					IJ.showStatus("distributing " + (fr+1) + "/" + impDepth);
				}
			}
			Odata.setStack(OdataSt);
			Fdata.setStack(FdataSt);
			for(int sfb = 0; sfb < sfbins.length; ++sfb)
			{
				Sfbin sfbin = sfbins[sfb];
				sfbin.imp.setStack(sfbin.impSt);
			}
			return true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return false;
		}
	}
	
	private boolean fetch_image(String title)
	{
		//System.out.println("ftitle: " + title + "\t imp: " + ((tmpImp==null)?"null":tmpImp.getTitle()));
		
		if (title!=null)
		{	
			tmpImp = WindowManager.getImage(title);
			if(tmpImp!=null)
			{	tmpImp.show();}	
		}
		
		if( (title!=null) && (tmpImp == null))
		{ 
			if(  ( !title.equals("<open>") )
			&& ( WindowManager.getImage(title) == null )  )
			{
				//IJ.log("could not open " + tmp.getPath());
				File tmp = new File(sourcepath + title + ".zip");
				if(tmp.isFile())
				{
					tmpImp = new ImagePlus (tmp.getPath());
					tmpImp.show();
				}
				else
				{
					tmp = new File(sourcepath + title);
					if(tmp.isFile())
					{
						tmpImp = new ImagePlus (tmp.getPath());
						tmpImp.show();
					}
					else
					{IJ.log("could not open " + sourcepath + title);}
				}
			}	
			
		}
		boolean status = ( (title!=null) && (tmpImp != null) && (WindowManager.getImage(title) == tmpImp )  );
		//System.out.println("status: " + status + "\timp: " + ((tmpImp==null)?"null":tmpImp.getTitle()));
		
		return status;
	}
	
	
	
	private boolean loadOptions()
	{
		try
		{
			if(options.isFile() && options.setReadable(true) )
			{
				IJ.log("reading " + options.getPath());
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
							{	Datatitle = ls.nextLine().trim();}
							else if(token.equals("#frame_name"))
							{	frame_name = ls.nextLine().trim();}
							else if(token.equals("#filter_name"))
							{	filter_name = ls.nextLine().trim();}
							else if(token.equals("#new_avg_level"))
							{	new_avg_level = ls.nextInt();}
							else if(token.equals("#divide_raw"))
							{	 divide_raw = ls.nextDouble();}
							else if(token.equals("#new_max_level"))
							{	new_max_level = ls.nextInt();}
							else if(token.equals("#min_histwidth"))
							{	min_histwidth = ls.nextInt();}
							else if(token.equals("#min_contrast"))
							{	min_contrast = ls.nextDouble();}
							else if(token.equals("#keep_Idata"))
							{	keep_Idata = ls.nextBoolean();}
							else if(token.equals("#prefix"))
							{	prefix = ls.nextLine().trim();}
							else if(token.equals("#bin_name"))
							{	bin_name = ls.nextLine().trim();}
							else if(token.equals("#unique_binning"))
							{	unique_binning = ls.nextBoolean();}
							else if(token.equals("#check_duplicates"))
							{	check_duplicates = ls.nextBoolean();}
							else if(token.equals("#Other_Data"))
							{	Othertitle = ls.nextLine().trim();}
							else if(token.equals("#do_dialog"))
							{	do_dialog = ls.nextBoolean();}
							else if(token.equals("#bondlength"))
							{	bondlength = ls.nextInt();}
							else if(token.equals("#modelGraylevels"))
							{	modelGraylevels = ls.nextInt();}
							else if(token.equals("#minfactor"))
							{	minfactor = ls.nextInt();}
							else if(token.equals("#modelBrightness"))
							{	modelBrightness = ls.nextDouble();}
							else if(token.equals("#gb_sigma"))
							{	gb_sigma = ls.nextDouble();}
							else if(token.equals("#pdf"))
							{	PDF = ls.nextLine().trim();}
							else if(token.equals("#pvalfloor"))
							{	pvalfloor = ls.nextDouble();}
							else if(token.equals("#interpolate"))
							{	interpolate = ls.nextBoolean();}
							else if(token.equals("#keep_step"))
							{	keep_step = ls.nextBoolean();}
							else if(token.equals("#sym_Histo"))
							{	sym_Histo = ls.nextBoolean();}
							else
							{	System.out.println("WARNING unknown token in options.txt: " + token + "\t" + ls.next());}
						}
					}
				}
			}
			
			return true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
			IJ.log("Error reading " + options.getPath() );
			return false;
		}
	}
	
	private boolean saveOptions(boolean tipps)
	{
		try
		{
			java.io.FileOutputStream opts = new java.io.FileOutputStream(options, false);
			java.io.PrintWriter optp = new java.io.PrintWriter(opts);
			optp.println("options are read from lines that start with #");
			optp.println("Please make sure your filenames do not contain any blank characters");
			optp.println("#appname\t" + getClass().getSimpleName() );
			optp.println("#sourcepath\t" + sourcepath);
			optp.println("#Input_Data\t" + Datatitle);
			optp.println("#frame_name\t" + frame_file.getName());
			optp.println("#filter_name\t" + filter_name);
			optp.println("#new_avg_level\t" + new_avg_level);
			optp.println("#divide_raw\t" + divide_raw);
			optp.println("#new_max_level\t" + new_max_level);
			optp.println("#min_contrast\t" + min_contrast);
			optp.println("#min_histwidth\t" + min_histwidth);
			optp.println("#keep_Idata\t" + keep_Idata);
			optp.println("#prefix\t" + prefix);
			optp.println("#bin_name\t" + bin_file.getName());
			optp.println("#unique_binning\t" + unique_binning);
			optp.println("#check_duplicates\t" + check_duplicates);
			if(tipps) optp.println("left over frames are stored here (will be over written)");
			optp.println("#Other_Data\t" + Othertitle);
			if(tipps) optp.println("if the main dialog should be shown, file openers are not suppressed");
			optp.println("#do_dialog\t" + do_dialog);
			optp.println("#bondlength\t" + bondlength);
			optp.println("#modelGraylevels\t" + modelGraylevels);
			optp.println("#modelBrightness\t" + modelBrightness);
			optp.println("#pdf\t" + PDF);
			optp.println("#minfactor\t" + minfactor);
			optp.println("#gb_sigma\t" + gb_sigma);
			optp.println("#pvalfloor\t" + pvalfloor);
			optp.println("#interpolate\t" + interpolate);
			optp.println("#keep_step\t" + keep_step);
			optp.println("#sym_Histo\t" + sym_Histo);
			
			
			
			
			
			optp.flush();
			optp.close();
			opts.close();
			return true;
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return false;
		}
	}
}
