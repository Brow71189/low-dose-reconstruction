import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.*;
import java.util.Arrays;
import java.awt.Point;

public class Frame_Aligner implements PlugIn 
{
	Point loc = null;
	
	//coords of active Hex
	int cq = 0;
	int cr = 0;
	int cx = cq - (cr - (cr&1)) / 2;
	int cz = cr;
	int cy = -cx - cz;
	int cradius = 0;//only the hex at the center
	int imgedge = 0;
	
	int mapW = -1;
	int mapH = -1;
	
	int search_radius = 4;
	int stepsize = 0;
	
	int firstframe = 1;
	int lastframe = -1;
	int latticeframe = 1;
	int bgframe = 1;
	
	int frames_done = 0;
	int moved_frames = 0;
	
	//effective units for faked Hex_Symmetry
	int hpMax = 1;
	int lpMax = 1;
	int hp000 = 0;
	int lp000 = 0;
	int transMax = 1;
	int mirrorMax = 1;
	int rotMax = 1;
	
	int[] hexPixels = null;
	int[] transX = null;
	int[] transZ = null;
	int[][] transformations = null;
	
	//these are used for extracting matches
	double min_corr = -1.0;
	double max_corr = 1.0;
	
	boolean align_slices = false;
	boolean full_hexagons = false;
	boolean update = true;
	boolean virtual_smoothing = false;
	boolean enable_rotations = false;
	boolean enable_mirror = false;
	boolean filter_matches = false;
	boolean enforce_multithreading = true;
	boolean ignore_zeros = false;
	
	//NOTHING	0;
	//DONE 4096;
	int flags = 0; 
	
	public Frame_Aligner(){}
	private OutlierCore outlierCore = null;
	private float[] matches = null;
	
	
	//normally all three arguments would be width/2
	public void setHexcenter(int q, int r, int rH)
	{
		//in principle negative centers are allowed
		cq = q;
		cr = r;
		cradius = rH;
		cx = q - (r-(r&1)) / 2;
		cz = r;
		cy = -cx - cz;
	}
	
	//performs close to center test for "odd-r" hexes 
	private boolean iswithinHex(int q, int r)
	{
		if(q < 0 || r < 0)
		{	IJ.error("Frame_Aligner.iswithinHex: image coordinates cannot be negative");}
		//cube coords of target hex
		int px = q - (r - (r&1)) / 2;
		int pz = r;
		int py = -px - pz;
		//offset hex in cube coords
		int dx = px - cx;
		int dy = py - cy;
		int dz = pz - cz;
		if(full_hexagons)
		{
			return ( (Math.abs(dx) + Math.abs(dy) + Math.abs(dz) ) <= 2 * cradius);	
		}
		else
		{
			return (dx >= -cradius && dx <  cradius && 
					dy >  -cradius && dy <= cradius &&
					dz >= -cradius && dz <  cradius);
		}	
	}
	
	//apply periodic boundary conditions to the RELATIVE coordinates
	public void hexagonaltrap(int[] cube, int N)
	{
		int lx = cube[0];
		int ly = cube[1];
		int lz = cube[2];
		//System.out.print( "x, y, z: " + lx + ", " + ly + ", " + lz);
		int alx = (int)Math.abs(lx);
		int aly = (int)Math.abs(ly);
		int alz = (int)Math.abs(lz);
							
		if(full_hexagons)
		{
			boolean is_inside = true;
			do
			{
				is_inside = true;
				//positive borders
				while (lx > N)
				{
					is_inside = false;
					lx -= 2*N + 1 ;
					ly += N;
					lz += N + 1;
				}
				while (ly > N)
				{
					is_inside = false;
					lx += N + 1;
					ly -= 2*N + 1;
					lz += N;
				}
				while (lz > N)
				{
					is_inside = false;
					lx += N;
					ly += N + 1;
					lz -= 2*N + 1;
				}
				//negative borders
				while (lx < -N)
				{
					is_inside = false;
					lx += 2*N + 1;
					ly -= N;
					lz -= N + 1;
				}
				while (ly < -N)
				{
					is_inside = false;
					lx -= N + 1;
					ly += 2*N + 1;
					lz -= N;
				}
				while (lz < -N)
				{
					is_inside = false;
					lx -= N;
					ly -= N + 1;
					lz += 2*N + 1;
				}	
				//System.out.println( " -> " + lx + ", " + ly + ", " + lz);
			}
			while(!is_inside);
		}
		else 
		{
			
			boolean at_home = true;
			//apply periodic boundary conditions
			do
			{
				at_home = true;
				if(lx >= N)
				{
					at_home = false;
					lx -= (2*N); //+1
					ly -= (-N);
					lz -= (-N); //-1
				}
				else if (lx < -N)
				{
					at_home = false;
					lx += (2*N); //+1
					ly += (-N);
					lz += (-N); //-1
				}
			
				if(ly > N)
				{
					at_home = false;
					lx -= (-N); //-1
					ly -= (2*N); //+1
					lz -= (-N);
				}
				else if (ly <= -N)
				{
					at_home = false;
					lx += (-N); //-1
					ly += (2*N); //+1
					lz += (-N);
				}
			
				if(lz >= N)
				{
					lx -= (-N);
					ly -= (-N); //-1
					lz -= (2*N); //+1
				}
				else if (lz < -N)// if(lz < 0)
				{
					lx += (-N);
					ly += (-N); //-1
					lz += (2*N); //+1
				}	
			} while (!at_home);
			
		}
		 
		cube[0] = lx;
		cube[1] = ly;
		cube[2] = lz;
	}
	
	private double hexAvg( short[] pixels )
	{
		double Sum = 0.0;
		final int rH = cradius;
		int non_zeros = 0;
		if(!ignore_zeros)
		{
			for(int hpC = 0; hpC < hpMax; ++hpC)
			{
				Sum += (double)pixels[hexPixels[hpC]];		
			}
			non_zeros = hpMax;
		}
		else
		{
			for(int ind = 0; ind < pixels.length; ++ind)
			{
				Sum += (double)pixels[ind];
				if(pixels[ind] != 0)
				{	++non_zeros;}		
			}
		
		}	
		return Sum/non_zeros;
	}
	
	private double hexAvg( double[] pixels )
	{
		double Sum = 0.0;
		final int rH = cradius;
		int non_zeros = 0;
		if(!ignore_zeros)
		{
			for(int hpC = 0; hpC < hpMax; ++hpC)
			{
				Sum += pixels[hexPixels[hpC]];		
			}
			non_zeros = hpMax;
		}
		else
		{
			for(int ind = 0; ind < pixels.length; ++ind)
			{
				Sum += pixels[ind];
				if(pixels[ind] != 0.0)
				{	++non_zeros;}		
			}
		
		}	
		return Sum/non_zeros;
	}
	
	
	
	private double hexSigma( short[] pixels, double mean )
	{
		double Sum2 = 0.0;
		int non_zeros = 0;
		if(!ignore_zeros)
		{
			for(int hpC = 0; hpC < hpMax; ++hpC)
			{
				final double s = pixels[hexPixels[hpC]]-mean;
				Sum2 += s*s;		
			}
			non_zeros = hpMax;
		}
		else
		{
			for(int ind = 0; ind < pixels.length; ++ind)
			{
				final double s = pixels[ind]-mean;
				Sum2 += s*s;
				if(pixels[ind] != 0)
				{	++non_zeros;}		
			}	
		}
		return Math.sqrt(Sum2/non_zeros);

	}
	
	private double hexSigma( double[] pixels, double mean )
	{
		double Sum2 = 0.0;
		int non_zeros = 0;
		if(!ignore_zeros)
		{
			for(int hpC = 0; hpC < hpMax; ++hpC)
			{
				final double s = pixels[hexPixels[hpC]]-mean;
				Sum2 += s*s;		
			}
			non_zeros = hpMax;
		}
		else
		{
			for(int ind = 0; ind < pixels.length; ++ind)
			{
				final double s = pixels[ind]-mean;
				Sum2 += s*s;
				if(pixels[ind] != 0)
				{	++non_zeros;}		
			}	
		}
		return Math.sqrt(Sum2/non_zeros);

	}
	

	private double hex_crossc(double[] lattice,  short[] source, double avg_source, int lp) //double avg_lattice,
	{
		double ProdSum = 0.0;
		int non_zeros = 0;
		if(!ignore_zeros)
		{
			for(int hp = 0; hp < hpMax; ++hp)
			{		
				final int indl = hexPixels[hp];
				final int inds = transformations[hp][lp];		
				ProdSum += ( ( lattice[indl]  ) * //- avg_lattice
							 ( (double)source[inds] - avg_source  ) ); 
			}
			non_zeros = hpMax;
		}
		else
		{
			for(int hp = 0; hp < hpMax; ++hp)
			{
				final int indl = hexPixels[hp];
				final int inds = transformations[hp][lp];
				ProdSum += ( ( lattice[indl] ) * // - avg_lattice
							 ( (double)source[inds] - avg_source  ) ); 
				if( (lattice[indl] != 0.0) || (source[inds]!=0) )
				{	++non_zeros;}		
			}	
		}
		
			
		return (ProdSum/non_zeros);
	}
	

	private class AlignThread extends Thread
	{
		int start;
		int end;
		int _frames_done = 0;
		int _moved_frames = 0;
		ImagePlus latticeimg;
		ImagePlus bgimg;
		ImagePlus sources[];
		
		public AlignThread(int start, int end, ImagePlus latticeimg, ImagePlus bgimg, ImagePlus sources[])
		{
			this.start = start;
			this.end = end;
			this.latticeimg = latticeimg;
			this.bgimg = bgimg;
			this.sources = sources;
		}
		
		public void run()
		{
			alignFrames(start, end, latticeimg, bgimg, sources, this);
		}
	}
	
	private boolean alignFrames(int start, int end, ImagePlus latticeimg, ImagePlus bgimg, ImagePlus sources[], AlignThread my_thread)
	{
		try
		{
			final boolean hasbg = (bgimg!=null);
			ImageStack sourceSt = sources[0].getStack();
			ImageStack latticeSt = latticeimg.getStack();
			ImageStack bgSt = hasbg?bgimg.getStack():null;
			short[] lattice = (short[])latticeSt.getPixels(latticeframe);
			double avg_lattice = hexAvg(lattice);
			double sig_lattice = hexSigma(lattice, avg_lattice);
			double avg_k = avg_lattice;
			double sig_k = sig_lattice;
			double[] k = new double[lattice.length];
			
			short[] bg = hasbg?(short[])bgSt.getPixels(bgframe):null;
			double avg_bg =hasbg?hexAvg(bg):Double.NaN;
			double sig_bg =hasbg?hexSigma(bg,avg_bg):Double.NaN; 
			short[] source = null;
			double cov_bmp1 = 0.0;
			if(hasbg)
			{
				final double isig = 1.0/(sig_lattice*sig_bg);
				final double isig_bg = 1.0/(sig_bg);
				double[] bgd = new double[k.length]; 
				for(int i = 0; i < k.length; ++i)
				{
					k[i] = (	sig_bg*((double)lattice[i]-avg_lattice) - 
								sig_lattice*((double)bg[i]-avg_bg)	      
						   ) * isig;
					bgd[i] = ((double)bg[i]-avg_bg)*isig_bg;	   
				}
				cov_bmp1 = 1 - hex_crossc(bgd, lattice, avg_lattice, lp000) / sig_lattice;
				avg_k = hexAvg(k);
				sig_k = hexSigma(k,avg_k);
			}
			
			//normalize k, it is either the lattice or the scaled
			//difference between model and background
			{
				final double isig_k = 1.0/sig_k;
				for(int i = 0; i < k.length; ++i)
				{
					k[i] = ((hasbg?k[i]:(double)lattice[i])-avg_k)*isig_k;
				}
				avg_k = 0.0;
				sig_k = 1.0;
			}
			
			//that should help JIT to optimize
			final double avg_kc = avg_k;
			final double sig_kc = sig_k;
			final double cov_bmp1c = cov_bmp1;
			
			for(int ss = start; ss <= end; ++ss )
			{
				source = (short[])sourceSt.getPixels(ss);
				if(virtual_smoothing)
				{
					source = Arrays.copyOf(source,source.length);
					outlierCore.run(source);
				}
				
				final double avg_source = hexAvg(source);
				final double sig_source = hexSigma(source, avg_source);
				final double isig = 1.0 / (sig_source); //sig_kc *
				double match = -1.0; //will overwrite with first value in ANY case
				double minor_match = -1.0; //2nd best to measure significance
				//robust brute force search of highest crosscorelation
				boolean first_val = true;
				int top_lp = lp000;
				for(int lp=0; lp<lpMax ; ++lp )
				{			
					final double cross = hex_crossc(k, source, avg_source, lp); //lattice, avg_kc, source
					final double crosscorr = cross * isig + cov_bmp1c;
					/*
					if( (cross * isig > 1.0001) || (cross * isig < -1.0001) )
					{
						throw new IllegalStateException("crosscorrelation out of range");
					}
					*/
					if( (first_val) || (crosscorr > match) ) //>= fails with perfect test data
					{
						top_lp = lp;
						minor_match = match;
						match = crosscorr;
						if(first_val)
						{
							first_val = false;
							minor_match = match;
						}
					}
				}
				matches[ss-1]=(float)match;
				if( !Double.isNaN(match) && (top_lp != lp000) ) //as good as && (match==match)
				{
					if(align_slices)
					{
						for(int i=0;i<sources.length; ++i) //and the mapped stacks as well
						{
							source = (short[])sources[i].getStack().getPixels(ss);
							short[] backup = Arrays.copyOf(source,source.length);
							for(int hp = 0; hp < hpMax; ++hp)
							{
								int inds = hexPixels[hp];
								int indb = transformations[hp][top_lp]; 
								source[inds] = backup[indb];
							}
						}
					}
					if(my_thread != null)
					{
						++my_thread._moved_frames; 
					}
					else
					{
						++moved_frames;
					}
				}
				if(my_thread != null)
				{
					++my_thread._frames_done; 
				}
				else
				{
					++frames_done;
				}
			}
			
			return true;
		} 
		catch (Exception e)
		{
			e.printStackTrace();
            return false;
		}
	}
	
	public void run(String arg)
    {
		boolean success = true;
		String sourcetitle = null;
		String valtitle1 = null;
		String valtitle2 = null;
		String valtitle3 = null;
		String bgtitle = null;
		String latticetitle = null;
		
		
		do //one time loop unless an input error causes continue
		{
			boolean bad_input = false;
			NonBlockingGenericDialog gd = new NonBlockingGenericDialog(getClass().getSimpleName());
			int[] idArray = WindowManager.getIDList(); // list of all opened images (IDs)	
			if( (idArray == null) ) //just to be on the safe side
			{
				IJ.error("no/too few open images, quitting " + getClass().getSimpleName());
				return;
			}
			String[] titleArray = new String[idArray.length]; // titles of opened images
			String[] titleArrayP = new String[idArray.length+1];
			for (int i = 0; i < idArray.length; ++i)
			{
				titleArray[i] = WindowManager.getImage(idArray[i]).getTitle();
				titleArrayP[i+1] = titleArray[i];
			}
			
			String candidate1 = titleArray[0]; //the oldest
			String candidate2 = "<none>";
			titleArrayP[0] = candidate2;
			
			if(sourcetitle == null) sourcetitle = candidate1;
			if(valtitle1 == null) valtitle1 = candidate2;
			if(valtitle2 == null) valtitle1 = candidate2;
			if(valtitle3 == null) valtitle1 = candidate2;
			if(latticetitle == null) latticetitle = candidate1;
			if(bgtitle == null) bgtitle = candidate2;
			
			gd.addMessage("Images must be hex sampled");
			gd.addCheckbox("virtual smoothing", virtual_smoothing);
			gd.addChoice("key Stack (GRAY16)", titleArray, sourcetitle );
			gd.addChoice("value Stack 1 (GRAY16)", titleArrayP, valtitle1 );
			gd.addChoice("value Stack 2 (GRAY16)", titleArrayP, valtitle2 );
			gd.addChoice("value Stack 3 (GRAY16)", titleArrayP, valtitle3 );
			gd.addNumericField("First frame to process", firstframe, 0);
			gd.addNumericField("Last frame to process (<1 End)", lastframe, 0);
			gd.addChoice("Ideal lattice or model (GRAY16)", titleArray, latticetitle );
			gd.addNumericField("Lattice/model slice", latticeframe, 0);
			gd.addChoice("common background(will disable aligning)", titleArrayP,bgtitle);
			gd.addNumericField("background slice", bgframe, 0);
			gd.addNumericField("searchradius)", search_radius, 0);
			gd.addNumericField("stepsize",stepsize,0);
			gd.addCheckbox("align slices", align_slices);
			gd.addCheckbox("rotations", enable_rotations);
			gd.addCheckbox("mirror", enable_mirror);
			gd.addCheckbox("filter matches", filter_matches);
			gd.addCheckbox("enforce multithreading", enforce_multithreading);
			gd.addCheckbox("ignore zeros or hexdomain", ignore_zeros);
			
			if(loc != null)
			{
				gd.centerDialog(false);
				gd.setLocation(loc.x+10,loc.y+10); //compensate the drift
			}
			gd.showDialog();
			if (gd.wasCanceled())
			{
				break;
			}
			loc = gd.getLocation();
			
			virtual_smoothing = gd.getNextBoolean();
			sourcetitle = gd.getNextChoice();
			valtitle1 = gd.getNextChoice();
			valtitle2 = gd.getNextChoice();
			valtitle3 = gd.getNextChoice();
			firstframe = (int)gd.getNextNumber();
			lastframe = (int)gd.getNextNumber();
			latticetitle = gd.getNextChoice();
			latticeframe = (int)gd.getNextNumber();
			bgtitle = gd.getNextChoice();
			bgframe = (int)gd.getNextNumber();
			search_radius = (int)gd.getNextNumber();
			stepsize = (int)gd.getNextNumber();
			align_slices = gd.getNextBoolean();
			enable_rotations = gd.getNextBoolean();
			enable_mirror = gd.getNextBoolean();
			filter_matches = gd.getNextBoolean();
			enforce_multithreading = gd.getNextBoolean();
			ignore_zeros = gd.getNextBoolean();
			
			ImagePlus sourceimg = WindowManager.getImage(sourcetitle);
			if(sourceimg == null)
			{
				IJ.error("Could not find " + sourcetitle);
				bad_input = true;
			}
			
			ImagePlus valimg1 = WindowManager.getImage(valtitle1);
			ImagePlus valimg2 = WindowManager.getImage(valtitle2);
			ImagePlus valimg3 = WindowManager.getImage(valtitle3);
			ImagePlus bgimg = WindowManager.getImage(bgtitle);
			
			ImagePlus latticeimg = WindowManager.getImage(latticetitle);
			if(latticeimg == null)
			{
				IJ.error("Could not find " + latticetitle);
				bad_input = true;
			}
			if( 	(stepsize > search_radius) ||
					(stepsize < 0) ||
					( (stepsize > 0) && (search_radius%stepsize != 0) )  )
			{
				IJ.error("searchradius and stepsize dont fit");
				bad_input = true;
			}
			
			if(bad_input)
			{	
				success = false;
				continue;
			} //really too bad input here
			int num_slices = sourceimg.getStackSize();
			
			mapW = (int)Math.sqrt(num_slices);
			mapH = (num_slices+mapW-1)/mapW;
			matches = new float[mapW*mapH];
			Arrays.fill(matches,Float.NaN);
			
			ImagePlus sources_raw[] = {sourceimg,valimg1,valimg2,valimg3};
			int num_sources = 0;
			for(int i = 0; i < sources_raw.length; ++i)
			{
				if(sources_raw[i] != null)
				{
					++num_sources;
					sources_raw[i].show();
				}
			}
			ImagePlus sources[] = new ImagePlus[num_sources];
			int k = 0;
			for(int i = 0; i < sources_raw.length; ++i)
			{
				if(sources_raw[i] != null)
				{
					sources[k++]=sources_raw[i];
				}
			}
			
			imgedge = latticeimg.getWidth();
			
			if(  imgedge != latticeimg.getHeight() )
			{
				IJ.error("Lattice/model must be a square image");
				bad_input = true;
			}
			
			if ( (latticeimg.getType() != ImagePlus.GRAY16))
			{
				IJ.error(latticetitle + " is not GRAY16 formated");
				bad_input = true;
			}
			
			if (bgimg!=null) 
			{
				if(bgimg.getType() != ImagePlus.GRAY16)
				{
					IJ.error(bgtitle + " is not GRAY16 formated");
					bad_input = true;
				}
				if( (bgimg.getWidth() != imgedge) || (bgimg.getHeight() != imgedge ) )
				{
					IJ.error(bgtitle + " has wrong dimensions");
					bad_input = true;
				}
			
			}
			for(int i = 0; i < num_sources; ++i) //no need to check the sourceimg again
			{
				if( (imgedge != sources[i].getWidth() ) || 
					(imgedge != sources[i].getHeight() ) ||
					(num_slices != sources[i].getStackSize() ) ||
					(sources[i].getType() != ImagePlus.GRAY16) )
				{
					IJ.error("wrong dimensions and/or format of " + sources[i].getTitle() );	
					bad_input = true;
				}
				for(int m=0; m < i; ++m) //make sure there are no duplicates
				{
					if(sources[m]==sources[i])
					{	
						IJ.log("Duplicate image found " + sources[m].getTitle());
						bad_input = true;
					}
					
				}
			}
			
			
			
			if(bad_input){	continue;}			
			if(virtual_smoothing)
			{
				if(outlierCore == null)
				{	outlierCore = new OutlierCore(flags);}
				outlierCore.create_report = true;
				virtual_smoothing = (outlierCore.showDialog(sources[0],"virtual Outlier Masking") == flags);
			}
			
				
			if(firstframe > lastframe && lastframe > 1)
			{
				int tmp = lastframe;
				lastframe = firstframe;
				firstframe = tmp;
			}
			
			if(firstframe < 1)
			{
				firstframe = 1;
			}
			
			if(lastframe > num_slices || lastframe < 1)
			{
				lastframe = num_slices;
			}
			
			if(latticeframe < 1)
			{	
				latticeframe = 1;
			}
			
			if(latticeframe > latticeimg.getStackSize())
			{
				latticeframe = latticeimg.getStackSize();
			}
			int threads = Runtime.getRuntime().availableProcessors();
			int num_frames = lastframe - firstframe + 1;
			if( (!enforce_multithreading) && (num_frames < 1000))
			{	threads = 1;}
			
			IJ.log(getClass().getSimpleName());
			IJ.log("Source: " + sourcetitle + "   slices " + firstframe  +" to " + lastframe); 
			IJ.log("Lattice: " + latticetitle + " slice: " + latticeframe);
			if(bgimg!=null)
			{
				IJ.log("Background: " + bgtitle + " slice: " + bgframe + (align_slices?" refusing actual alignment":"") );
				align_slices = false;
			}	 
			IJ.log("search_radius: " + search_radius + " step_size_ " + stepsize +  " threads: " + threads);
			
			if(enable_mirror || enable_rotations || align_slices )
			{	IJ.log("mirror: " + enable_mirror + "   rotations: " + enable_rotations + " align slices: " + align_slices);}
			
			long timeoffset = System.currentTimeMillis();
			frames_done = 0;
			moved_frames = 0;
			
			setHexcenter(imgedge/2, imgedge/2, imgedge/2);
			init_hexPixels(); //provides lookup looping over frames
			init_transformations(); //provides lookup translations (also mirrors, rotations)
			
			if( (threads == 1) && (num_frames < 256) ) 
			{   //we dont need any progressbar or esc interrupt for quick jobs
				alignFrames(firstframe, lastframe, latticeimg, bgimg, sources, null);
			}
			else 
			{   //multi threading with interrups even for single thread if the job is big enough
				int assigned = 0;
				int chunk_start[] = new int[threads];
				int chunk_end[] = new int[threads];
				for(int i=0; i < threads; ++i)
				{
					int portion = (num_frames-assigned)/(threads-i);
					chunk_start[i] = assigned + 1; //stupid one based indexing
					chunk_end[i] = (assigned += portion); //requires inclusive test in for loop
				}

				AlignThread at[] = new AlignThread[threads]; 
				for (int i=0; i<threads; ++i)
				{
					//IJ.log("thread: " + i + " first: " + chunk_start[i] + " last: " + chunk_end[i]);
					at[i] = new AlignThread(chunk_start[i],chunk_end[i],latticeimg, bgimg, sources);
					at[i].start();
				}
				int joined = 0;
				try
				{
					do
					{
						int active = 0;
						frames_done = 0;
						moved_frames = 0;
						for(int i = 0; i < threads; ++i)
						{
							frames_done += at[i]._frames_done;
							moved_frames += at[i]._moved_frames;
							if(at[i].isAlive()) 
							{++active;}
						}
						IJ.showProgress(frames_done,num_frames);
						IJ.showStatus(getClass().getSimpleName() + " " + frames_done + "/" + num_frames);
						if(IJ.escapePressed())
						{
							IJ.log("Thread alignment was interrupted by Esc");
							IJ.resetEscape();
							break;
						}
						Thread.sleep(250);
						joined = threads - active;
					}
					while(joined < threads);
				}
				catch(Exception e)
				{
					e.printStackTrace();
					//success = false; //we should still be fine if the workes have completed their tasks
				}
				finally
				{
					frames_done = 0;
					moved_frames = 0;
					for(int i = 0; i < threads; ++i)
					{
						frames_done += at[i]._frames_done;
						moved_frames += at[i]._moved_frames;
						
						if(at[i].isAlive())
						{	
							at[i].interrupt();
							success = false;
						} //very brutal, but that is what we want here
					}
					IJ.showProgress(frames_done,num_frames);
					IJ.showStatus(getClass().getSimpleName() + " " + frames_done + "/" + num_frames);
				}
			} //end multi threading
			
			for(int i=0; i < sources.length; ++i)
			{	sources[i].updateAndRepaintWindow();}
			
			//sourceimg.updateAndRepaintWindow();
			
			IJ.log( (align_slices?"moved ":"missaligned slices: ") + moved_frames + " of " + frames_done );
			IJ.log("Elapsed time: " + ((System.currentTimeMillis()-timeoffset)/1000) + "s" );
			if(filter_matches)
			{	run_matches(sourceimg);}
			
			if(!success)
			{
				IJ.log("ERROR in alignSubFrames");
			}
				
			
		}while(true); // !success
	}
	
	private void run_matches(ImagePlus sourceimg)
	{
		ImagePlus matchmap = NewImage.createFloatImage( WindowManager.getUniqueName("match_map"), mapW , mapH , 1, NewImage.FILL_BLACK);
		ImageStack matchSt = matchmap.getStack();
		matchSt.setPixels(matches,1);
		matchmap.setStack(matchSt);
		matchmap.show();
		boolean keep_going = true;
		boolean inverse_test = false;

		do
		{
			NonBlockingGenericDialog gd = new NonBlockingGenericDialog(getClass().getSimpleName() + " extracting slices");
			gd.addMessage("active match map: " + matchmap.getTitle() );
			gd.addNumericField("minimal correlation: ", min_corr, 4 );
			gd.addNumericField("maximal correlation: ", max_corr, 4 );
			gd.addCheckbox("invert range", inverse_test);
			
			gd.addCheckbox("interactive mode: ", keep_going );
			if(loc != null)
			{
				gd.centerDialog(false);
				gd.setLocation(loc.x+10,loc.y+10); //compensate the drift
			}
			
			gd.showDialog();
			if (gd.wasCanceled())
			{
				keep_going = false;
				break;
			}
			loc = gd.getLocation();
			min_corr = gd.getNextNumber();
			max_corr = gd.getNextNumber();
			inverse_test = gd.getNextBoolean();
			
			keep_going = gd.getNextBoolean();
			
			if(max_corr >= min_corr)
			{
				extract_frames(sourceimg,inverse_test);
			}

		}while (keep_going);		
	}
	
	private void extract_frames(ImagePlus sourceimg,boolean inverse_test)
	{
		ImagePlus extd = NewImage.createShortImage( WindowManager.getUniqueName("extracted"), imgedge , imgedge , 1, NewImage.FILL_BLACK);
		ImageStack extdSt = extd.getStack();
		ImageStack sourceSt = sourceimg.getStack();
		int slc = 0;
		
		for(int ss = firstframe; ss <= lastframe; ++ss)
		{
			if(Float.isNaN(matches[ss-1]))
			{	continue;} //skip nans
			boolean in_range = (matches[ss-1] > min_corr) && (matches[ss-1] <= max_corr);
			
			if( inverse_test ^ in_range ) // as good as !=
			{
				String label = sourceSt.getSliceLabel(ss);
				short[] pixels = (short[])sourceSt.getPixels(ss);
				pixels = Arrays.copyOf(pixels, pixels.length);
				extdSt.addSlice(label, pixels);
				++slc;
				if( (slc == 2) && (extdSt.getSize() == 3) )
				{	extdSt.deleteSlice(1);}
			}
		}
		extd.setStack(extdSt);
		extd.show();
		IJ.log("extracted " + extdSt.getSize() + " matches " + (inverse_test?"outside ": "inside ") + min_corr + " <= corr <= " + max_corr + " to " + extd.getTitle());	
	}
	
	
	
	private void init_hexPixels()
	{
		hpMax = ( 3 * imgedge * imgedge) / 4; // pixels per hex domain
		hexPixels = new int[hpMax];
		mirrorMax = enable_mirror?2:1;
		rotMax = enable_rotations?6:1;
		transMax = 3*search_radius*search_radius;
		if(stepsize > 0)
		{	transMax = ( search_radius * search_radius) / ( stepsize * stepsize ); }
		lpMax = mirrorMax*rotMax*transMax;
		transX = new int[transMax];
		transZ = new int[transMax];
		
		final int rH = cradius;
		int hpC = 0;
		int trC = 0;
		for(int dx = -rH; dx < rH ; ++dx)
		{
			int dz1 = -rH;
			int dz2 = rH - dx;
			if (dx < 0)
			{
				dz1 = -rH - dx;
				dz2 = rH;
			}
			for(int dz = dz1; dz < dz2; ++dz)
			{
				int dy = -dx - dz;
				//translate to center
				int x_data = dx + cx;
				int z_data = dz + cz;
				// cube -> odd-r
				int q_data = x_data + ( z_data-(z_data&1) ) / 2;
				int r_data = z_data;
				int ind = q_data + r_data * imgedge;
				hexPixels[hpC++] = ind;
				
				
				if(dx >= -search_radius && dx <  search_radius && 
					dy >  -search_radius && dy <= search_radius &&
					dz >= -search_radius && dz <  search_radius)
				{ //Yes we are inside the search radius
					
					boolean trans_ok = (stepsize==0); //any translation is valid with trivial stepsize
					if(!trans_ok) //lets check if the translation matches stepsize 
					{
						int cube[] = {dx,dy,dz};
						hexagonaltrap(cube,stepsize);
						trans_ok = (cube[0]==0 && cube[1]==0);
					}
						
					if(trans_ok)
					{
						transX[trC]=dx;
						transZ[trC]=dz;
						if( (dx==0) && (dz==0) )
						{	hp000 = trC;}
						trC++;	
					}
	
				}		
			}
		}
		//IJ.log("precached " + hpC + " pixels inside hexagonal domain and " + trC + " translations inside search radius");
		//IJ.log("hp000: " + hp000);	
	}
		
	private void init_transformations()
	{
		transformations = new int[hpMax][lpMax];
		int transC = 0;
		final int rH = cradius;
		for(int hpC=0; hpC < hpMax; ++hpC)
		{
			int lpC = 0;
			int ind0 = hexPixels[hpC];
			int q0 = ind0 % imgedge;
			int r0 = ind0 / imgedge;
			//offset around center
			int x0 = q0 - (r0-(r0&1)) / 2 - cx;
			int z0 = r0 - cz;
			int y0 = (-x0 - z0);
			
			for(int mirror = 0; mirror < mirrorMax; ++mirror)
			{
				//mirror on x-axis //actually flipping upside down in 3D
				if( (mirror&1) == 1)
				{
					z0 = y0;
					y0 = (-x0 - z0);
				}
				
				for(int rot = 0; rot < rotMax; ++rot)
				{
					//60° rotation
					if(rot > 0)
					{
						int[] rotxyz = {x0, y0, z0};
						x0 = rotxyz[rot % 3];
						y0 = rotxyz[(1+rot) % 3];
						if( (rot&1) == 1)
						{
							x0 *= -1;
							y0 *= -1;
						}
						z0 = (-x0 - y0);
					}
				
					for(int hpT=0; hpT < transMax; ++hpT)
					{
						if( (hpT == hp000) && (rot == 0) && (mirror==0) )
						{	lp000 = lpC;}
						//translate
						int xt = transX[hpT] + x0;
						int zt = transZ[hpT] + z0;
						int yt = -xt -zt;
						
						int[] cube = {xt,yt,zt};
						hexagonaltrap(cube,rH);
						//undo offset around center
						int xf = cube[0]+cx;
						int yf = cube[1]+cy;
						int zf = cube[2]+cz;
			
						// cube -> odd-r
						int qf = xf + (zf - (zf&1) ) / 2;
						int rf = zf;
						
						int indf = qf + rf * imgedge; 
						transformations[hpC][lpC++] = indf;
						++transC;
									
					} //end step z		
				} //end rot
			} //end mirror
		} //end hpC
		//IJ.log("precached " + transC + " transformations");
		//IJ.log("lp000: " + lp000);
	}

}
