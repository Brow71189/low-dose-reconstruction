import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.*;
import java.awt.Point;
//import java.util.StringBuilder;

public class Diff_Matcher implements PlugIn 
{

	//coords of active Hex
	int cq = 0;
	int cr = 0;
	int cx = cq - cr / 2;
	int cz = cr;
	int cy = -cx - cz;
	int cradius = 0;//only the hex at the center
	int imgedge = 0;
	int bondlength = 4; //needed in case Hex_Symmetry is employed
	int grayvalues = 0; //need for optional subtraction of median by Hex_Symmetry 
	
	int firstframe = 1;
	int lastframe = -1;
	
	boolean show_matches = false;
	boolean verboose = false;
	boolean update = false; //THIS MUST ALWAYS BE FALSE BEFORE run(String)!!!
	boolean previous = false;
	boolean use_symmetrizers = false; //THIS MUST ALWAYS BE FALSE BEFORE run(String)!!!
	boolean show_symImgs = false;
	boolean detached = false;
	boolean silent = false;
	public double summed_matches[] = null; //from here the resullts can be accessed
	public double inv_avg_matches[] = null;
	public double exponent = 1.0;
	
	String previous_best = "";
	int previous_best_pos = 1;
	String matchtitle = null;
	ImagePlus moleculeimg = null;
	ImagePlus diffimg = null;
	
	private Hex_Symmetry[] symmetrizers = null;
	
	public Diff_Matcher(){}
	
	//grayvalue > 2 will trigger median subtraction
	public Diff_Matcher(ImagePlus _diffimg, int _firstframe,
						int _lastframe, int _bondlength, int _grayvalues)
	{
		diffimg = _diffimg;
		imgedge = diffimg.getWidth();
		firstframe = _firstframe;
		lastframe = _lastframe;
		bondlength = _bondlength;
		grayvalues = _grayvalues;
		
		use_symmetrizers = true;
		update = true;
		show_symImgs = false; //we definitely dont want that
		verboose = false;
		detached = true; //overrides other settings anyways
		
		regulate_framerange();
		init_symmetrizers();
		update_matches(); //get the first summed_matches
		IJ.log("Instantiated a Diff_Matcher for " + diffimg.getShortTitle() + ( diffimg.isVisible()?"":"   (hidden)" ));
	}
	
	public void update_matches()
	{
		update_symmetrizers(); //update symmetrized stacks of the current diffimage
		matchMolecule(); //populate the public summed_matches
	}
	
	
	//normally all three arguments would be imgedge/2
	private void setHexcenter(int q, int r, int rH)
	{
		//in principle negative centers are allowed
		cq = q;
		cr = r;
		cradius = rH;
		cx = q - r / 2;
		cz = r;
		cy = -cx - cz;
	}
	
	//performs close to center test for "odd-r" hexes 
	private boolean iswithinHex(int q, int r)
	{
		if(q < 0 || q < 0)
		{	IJ.error("Diff_Matcher.iswithinHex: image coordinates cannot be negative");}
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
				dz >= -cradius && dz <  cradius);
			
	}
	
	//apply periodic boundary conditions to the RELATIVE coordinates
	private void hexagonaltrap(int[] cube, int N)
	{
		int lx = cube[0];
		int ly = cube[1];
		int lz = cube[2];
		//System.out.print( "x, y, z: " + lx + ", " + ly + ", " + lz);
		int alx = (int)Math.abs(lx);
		int aly = (int)Math.abs(ly);
		int alz = (int)Math.abs(lz);
								
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
			
		cube[0] = lx;
		cube[1] = ly;
		cube[2] = lz;
	}
	//TODO use "Verschiebungssatz" and populate a stadistics array
	private double hexAvg( short[] pixels )
	{
		int Npoints = 0;
		double Sum = 0.0;
		final int rH = cradius;
		for(int dx = -rH; dx <= rH ; ++dx) 
		{
			for(int dz = Math.max(-rH,-rH-dx); dz <= Math.min(rH, rH-dx); ++dz)
			{ 
				int dy = -dx - dz;
				if( (dz==rH) || (dy==-rH) || (dx==rH) )
				{
					continue;
				}
				int px = dx + cx;
				int pz = dz + cz;
				int q = (px + pz / 2);	
				int r = pz;
				
				++Npoints;
				Sum += (double) (pixels[q + r * imgedge] & 0xffff);
			}
		}
		return (Sum/Npoints);
	}
	
	private double hexSigma( short[] pixels, double mean )
	{
		int Npoints = 0;
		double Sum = 0.0;
		final int rH = cradius;
		for(int dx = -rH; dx <= rH ; ++dx) 
		{
			for(int dz = Math.max(-rH,-rH-dx); dz <= Math.min(rH, rH-dx); ++dz)
			{ 
				int dy = -dx - dz;
				if( (dz==rH) || (dy==-rH) || (dx==rH) )
				{
					continue;
				}
				int px = dx + cx;
				int pz = dz + cz;
				int q = (px + pz / 2);	
				int r = pz;
				
				++Npoints;
				Sum += Math.pow( (pixels[q + r * imgedge] & 0xffff ) - mean, 2);
			}
		}
		return Math.sqrt(Sum/Npoints);
	}
		
	private double hexcrosscorelate(short[] lattice, double avg_lattice, short[] source, double avg_source)
	{
		//int ty = -tx - tz;
		int Npoints = 0;
		double ProdSum = 0.0;
		final int rH = cradius;
		for(int dx = -rH; dx <= rH ; ++dx) 
		{
			for(int dz = Math.max(-rH,-rH-dx); dz <= Math.min(rH, rH-dx); ++dz)
			{ 
				int dy = -dx - dz;
				if( (dz==rH) || (dy==-rH) || (dx==rH) )
				{
					continue;
				}
				int px = dx + cx;
				int pz = dz + cz;
				int q = px + pz / 2;
				++Npoints;
				ProdSum += ( ( (double)(lattice[q + pz * imgedge] & 0xffff ) - avg_lattice ) * 
						     ( (double)(source[q + pz * imgedge] & 0xffff ) - avg_source  ) ); //int r = pz;
			}
		}
		return (ProdSum/Npoints);
	}
	
	
	private boolean matchMolecule()
	{
		try
		{
			previous_best = "Exception";
			int sym_n = 0;
			ImageStack diffSt = diffimg.getStack();
			ImageStack moleculeSt = use_symmetrizers ? symmetrizers[sym_n].OutputImg.getStack() : moleculeimg.getStack();
			int imgedge = diffimg.getWidth();
			ImageStack matchSt = show_matches ?
					ImageStack.create(imgedge,imgedge,diffimg.getStackSize(),16) :
					null;
			setHexcenter(imgedge/2, imgedge/2, imgedge/2);
			short[] molpix = (short[])moleculeSt.getPixels(1);
			
			double avg_molecule = hexAvg(molpix);
			double sig_molecule = hexSigma(molpix, avg_molecule);
			int mframes = moleculeSt.getSize();
			int[] best_pos = new int[diffSt.getSize()];
			double max_match = -1.0;
			
			int best_diff = 0;
			int worst_diff = 0;
			boolean median_filtered = false;
			boolean coherent = true;
			for(int ss = firstframe; ss <= lastframe; ++ss )
			{
				
				
				
				short[] diffpix = (short[])diffSt.getPixels(ss);
				if(symmetrizers != null && symmetrizers[ss-firstframe].median != null)
				{   //we need to subtract the median without changing the inputImg
					median_filtered = true;
					short[] ddiff = new short[diffpix.length];
					for(int ind = 0; ind < diffpix.length; ++ind)
					{
						ddiff[ind] = (short)(diffpix[ind] + 32768 - symmetrizers[ss-firstframe].median[ind]);
					}
					diffpix = ddiff; //median shifted private copy
				}
				else if (median_filtered)
				{
					coherent = false;
				}
				
				double avg_diff = hexAvg(diffpix);
				double sig_diff = hexSigma(diffpix, avg_diff);
				double sig_prod = 1.0 / (sig_molecule * sig_diff);
				boolean first_val = true;
				double match = -1.0; //will overwrite with first value in ANY case
				double anti_match = 1.0;
				double minor_match = -1.0; //2nd best to measure significance
				double sum_match = 0.0;
				for(int mfr = 1; mfr <= mframes; ++mfr)
				{	
					molpix = (short[])moleculeSt.getPixels(mfr);
					double crosscorr = hexcrosscorelate(molpix, avg_molecule, diffpix, avg_diff);
					
						crosscorr *= sig_prod;
						/**
						 * What does a negative correlation mean in terms  
						 * of symmetry breaking and likelyhoods? Should never
						 * happen anyways.
						
						**/
						if( crosscorr >= 0.0 ) 
						//also filters NaNs
						{	
							sum_match += crosscorr;
						}
						if( ( first_val && (crosscorr == crosscorr) ) || (crosscorr >= match) )
						{
							best_pos[ss-1] = mfr;
							minor_match = match;
							match = crosscorr;
							if(show_matches)
							{	matchSt.setPixels(molpix,ss);}
							if(first_val)
							{
								first_val = false;
								minor_match = match;
								anti_match = match;
							}
						}
						if( (crosscorr == crosscorr) && (crosscorr <= anti_match) )
						{
							anti_match = crosscorr;
						}
					
				}
				summed_matches[sym_n] = sum_match;
				inv_avg_matches[sym_n] = Math.pow(mframes/sum_match,exponent);
				if(!detached)
				{
					if(!update)
					{
						if(verboose)
						{	
							IJ.log("" + diffimg.getTitle() + " slice: " + ss + "  match: " + match + "  anti_match: " + anti_match);
							IJ.log("sum of (positive) matches: " + sum_match + "   avg. match: " + (sum_match/(mframes)));	
						}
						else if (!silent)
						{	IJ.log("" + match/* + "     " + anti_match*/);}
					}
				}
				if(ss == firstframe || match > max_match)
				{
					max_match = match;
					best_diff = ss;
				}
				if(use_symmetrizers && sym_n < (symmetrizers.length-1))
				{
					moleculeSt = symmetrizers[++sym_n].OutputImg.getStack();
					molpix = (short[])moleculeSt.getPixels(1);
					avg_molecule = hexAvg(molpix);
					sig_molecule = hexSigma(molpix, avg_molecule);
					mframes = moleculeSt.getSize();	
				}
			}
			if(!detached)
			{ 
				if(show_matches)
				{
					if(coherent && median_filtered)
					{
						for(int sl = 1; sl < matchSt.getSize(); ++sl)
						{
							short[] mpix = (short[]) matchSt.getPixels(sl);
							for(int px = 0; px < mpix.length; ++px)
							{
								int q = px % imgedge;
								int r = px / imgedge;
								if( iswithinHex(q, r) )
								{
									mpix[px] += 32768;
								}
							}
						}
					}
					
					
					ImagePlus matchimg = update ? WindowManager.getImage( (matchtitle==null) ? "matches_"+diffimg.getTitle() : matchtitle ) : null;
					if(matchimg == null)
					{
						matchimg = new ImagePlus(WindowManager.makeUniqueName("matches_"+diffimg.getTitle()),matchSt);
					}
					else
					{
						matchimg.setStack(matchSt);
					}
					matchimg.show();
					matchtitle = matchimg.getTitle();
				}
				if(verboose && !update)
				{
					IJ.log("best slice: " + best_diff + "  match: " + max_match);
				}
				if(update)
				{
					previous = true;
					if(!use_symmetrizers)
					{
						previous_best = "best slice: " + best_diff + "  match: " + max_match;
						previous_best_pos = best_pos[best_diff-1];
					}
					else
					{
						StringBuilder sb = new StringBuilder();
						if(coherent && !median_filtered)
						{
							sb.append("inv avg matches^" + exponent);
							for(int ss = firstframe; ss <= lastframe; ++ss)
							{
								sb.append("\nframe: " + ss + "  asym: " + inv_avg_matches[ss-firstframe] );
							}
						}
						else if(coherent && median_filtered)
						{
							sb.append("summed matches (max: " + mframes + ")");
							for(int ss = firstframe; ss <= lastframe; ++ss)
							{
								sb.append("\nframe: " + ss + "  sum: " + summed_matches[ss-firstframe] );
							}
						}
						else
						{	
							sb.append("incoherent median filtering!");
						}
						previous_best = sb.toString();
					}
				}
			}
			return true;
		} catch (Exception e)
		{
			e.printStackTrace();
            return false;
		}	
	}
	
	private void init_symmetrizers()
	{
		int num_frames = lastframe - firstframe + 1;
		symmetrizers = new Hex_Symmetry[num_frames];
		
		for(int nsym = 0; nsym < num_frames; ++nsym)
		{
			int fr = firstframe + nsym;
			symmetrizers[nsym] = 
				new Hex_Symmetry( fr, fr, imgedge, bondlength, grayvalues, diffimg);
		}
	} 
	
	private void update_symmetrizers()
	{
		int num_frames = lastframe - firstframe + 1;
		for(int nsym = 0; nsym < num_frames; ++nsym)
		{
			//update ImagePlus symmetrizers[nsym].OutputImg
			symmetrizers[nsym].apply_transformations();
			summed_matches[nsym] = Double.NaN; //they literally dont exist yet
			inv_avg_matches[nsym] = Double.NaN;
			if(show_symImgs)
			{
				symmetrizers[nsym].OutputImg.setTitle(diffimg.getShortTitle() + "_sym_" + (firstframe+nsym));
				symmetrizers[nsym].OutputImg.show();
			}
		}		
	}
	
	private void regulate_framerange()
	{
		int num_slices = diffimg.getStackSize();
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
		int num_frames = lastframe - firstframe + 1;
		summed_matches = new double[num_frames];
		inv_avg_matches = new double[num_frames];
	
	}
	
	
	public void run(String arg)
    {
		boolean success = true;
		String difftitle = "";
		String moleculetitle = "";
		Point loc = null;
		
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
			String[] titleArray2 = new String[idArray.length + 1];
			titleArray2[0] = "<symmetrize>";
			for (int i = 0; i < idArray.length; ++i)
			{
				titleArray[i] = WindowManager.getImage(idArray[i]).getTitle();
				titleArray2[i+1] = titleArray[i];
			}
			
			String candidate1 = "Diff.tif"; 
			
			String candidate2 = titleArray[idArray.length-1];//most recent
			
			if(difftitle.equals("") || (WindowManager.getImage(difftitle) == null) )
			{	difftitle = candidate1;}
			
			if(moleculetitle.equals("") || (WindowManager.getImage(moleculetitle) == null) )
			{	moleculetitle = titleArray2[0];}
			
			//gd.addMessage("Diff images must be hex sampled");
			if(!update)
			{	
				gd.addChoice("Diff Stack (GRAY16)", titleArray, difftitle );
				
				gd.addNumericField("First frame to process", firstframe, 0);
				gd.addNumericField("Last frame to process (<1 End)", lastframe, 0);
				//gd.addMessage("Stack of equivalent configurations");
				gd.addChoice("known molecule (GRAY16 stack)", titleArray2, moleculetitle );
				gd.addNumericField("bondlength",bondlength,0);
				gd.addNumericField("grayvalues(<2 no median filtering)", grayvalues, 0);
			}
			else
			{
				gd.addMessage("Diff Stack (GRAY16): " + difftitle);
				gd.addMessage("First Frame: " + firstframe);
				gd.addMessage("Last Frame: " + lastframe);
				gd.addMessage("known molecule (GRAY16): " + moleculetitle);
				gd.addMessage("bondlength: " + bondlength);
				gd.addMessage("gayvalues: " + grayvalues);
			}
			gd.addCheckbox("show matches", show_matches);
			gd.addCheckbox("no logging", silent);
			gd.addCheckbox("verboose logging", verboose);
			gd.addCheckbox("update", update);
			gd.addCheckbox("show symmetrized diffs", show_symImgs);
			
			if(previous)
			{
				gd.addMessage(previous_best);
				if(moleculeimg != null && previous_best_pos!=0)
				{	moleculeimg.setPosition(previous_best_pos);}
				previous = false; //only show this once
			}
			if(loc != null)
			{
				gd.centerDialog(false);
				gd.setLocation(loc.x+10,loc.y+10); //compensate the drift
			}
			gd.showDialog();
			loc = gd.getLocation(); 
			if (gd.wasCanceled())
			{
				break;
			}
			
			/*
			short a1 = (short)32766;
			short a2 = (short)32774;
			IJ.log("a1(32766): " + a1 + "  converted: " + ((int)(a1&0xffff)));
			IJ.log("a1(32774): " + a2 + "  converted: " + ((int)(a2&0xffff)));
			*/
			boolean fresh_input = false;
			if(!update)
			{
				fresh_input = true;
				difftitle = gd.getNextChoice();
				firstframe = (int)gd.getNextNumber();
				lastframe = (int)gd.getNextNumber();
				moleculetitle = gd.getNextChoice();
				bondlength = (int)gd.getNextNumber();
				grayvalues = (int)gd.getNextNumber();
			}
			show_matches = gd.getNextBoolean();
			
			silent = gd.getNextBoolean();
			verboose = gd.getNextBoolean();
			update = gd.getNextBoolean();
			show_symImgs = gd.getNextBoolean();
			
			if(verboose)
			{
				silent = false;
			}
			
			if(update && verboose)
			{
				IJ.showMessage("Hey dude!","verboose AND updating is no good idea");
				verboose = false;
			}
			if(fresh_input)
			{
				diffimg = WindowManager.getImage(difftitle);
				if(diffimg == null)
				{
					IJ.error("Could not find " + difftitle);
					bad_input = true;
				}
				if(!moleculetitle.equals("<symmetrize>"))
				{
					moleculeimg = WindowManager.getImage(moleculetitle);
					if(moleculeimg == null)
					{
						IJ.error("Could not find " + moleculetitle);
						bad_input = true;
					}
					imgedge = moleculeimg.getWidth();
				
					if(  imgedge != moleculeimg.getHeight() )
					{
						IJ.error("molecule must be a square image");
						bad_input = true;
					}
				
					if ( (moleculeimg.getType() != ImagePlus.GRAY16))
					{
						IJ.error(moleculetitle + " is not GRAY16 formated");
						bad_input = true;
					}	
				}
				else
				{	
					imgedge = diffimg.getWidth();
					use_symmetrizers = true;	
				}
				if(bad_input){	continue;} //really too bad input here
				
				
				
				if( (imgedge != diffimg.getWidth()  ) || 
					(imgedge != diffimg.getHeight() ) 	)  	
				{
					IJ.error("wrong dimensions and/or format of " + difftitle );	
					bad_input = true;
				}
				
				if(bad_input){	continue;}			
				
				regulate_framerange();
				if( use_symmetrizers )
				{
					init_symmetrizers();
				}
			
			}
			
			
			if(!update)
			{	
				IJ.log("Diff_Matcher\nDiffs: " + difftitle + "   frames " + firstframe  +" to " + lastframe + "\nMolecule: " 
					+ moleculetitle);
			}
			long timeoffset = System.currentTimeMillis();
			
			if( use_symmetrizers )
			{
				update_symmetrizers();
			}
			else if (show_symImgs)
			{
				IJ.log("There are no symmetrized Diffs to show");
				show_symImgs = false;
			}
			
			
			matchMolecule();
			
			if(verboose)
			{	IJ.log("Elapsed time: " + ((System.currentTimeMillis()-timeoffset)/1000) + "s" );}
			if(!success)
			{
				IJ.log("ERROR in matchMolecule");
			}
				
		}while(update || (!success) );
	}

}
