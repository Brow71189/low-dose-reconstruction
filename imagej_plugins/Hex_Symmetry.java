import ij.*;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.gui.GenericDialog;
import ij.plugin.*;
import java.util.Random;

public class Hex_Symmetry implements PlugIn {

	int Mx,Mz;
	int subframesize;
	int[] lpX = null;
	int[] lpZ = null;
	int[] hexPixels = null;
	int[] latticePoints = null;
	int hpMax = 1;
	int lpMax = 1; 
	int bondlength = 4;
	int subframeCount =1;
	int current_subframesize = -1;
	int current_bondlength = -1;
	int grayvalues = 64;
	
	String Inputtitle = WindowManager.getCurrentImage().getTitle(); //the lattice
	String Outputtitle = WindowManager.makeUniqueName("Hex_Symmetry.tif"); //the output
	public ImagePlus InputImg = null;
	public ImagePlus OutputImg = null;
	int firstframe = 1;
	int lastframe = -1;
	boolean quit_after_generation = false;
	boolean detached = false;
	boolean calculate_median = false;
	boolean show_median = false;
	boolean subtract_median = false;
	
	
	
	private short[][] full_histo = null;
	public short[] median = null; //needed by calling Diff_Matcher
	
	public Hex_Symmetry() {}
	
	public Hex_Symmetry(int _firstframe, int _lastframe, int _subframesize,
						int _bondlength, int _grayvalues, ImagePlus _InputImg)			 
	{
		detached = true;
		InputImg = _InputImg;
		firstframe = _firstframe;
		lastframe = _lastframe;
		subframesize = _subframesize;
		bondlength = _bondlength;
		grayvalues = _grayvalues;
		current_subframesize = _subframesize;
		current_bondlength = _bondlength;
		
		if(grayvalues > 1)
		{
			calculate_median = true;
			show_median = false;
			subtract_median = true;	
		}
		
		if(validateInput() < 2)
		{	init_transformations();}
		else
		{	System.out.println("ERROR in Hex_Symmetry::Hex_Symmetry(int, int, int, int, ImagePlus)"); }
		
		//IJ.log("new Hex_Symmetry instantiated");
		//IJ.log(InputImg.getTitle() + " first_frame: " + firstframe + "  last_frame: " + lastframe);
		//IJ.log("hidden " OutputImg.getTitle() + "  frames: " + OutputImg.getStackSize());
		
		/*	After constructing every call of apply_transformations()
		 *	will update OutputImg from the current state of the
		 *  input image
		 */  	
		
	}
	
	public void run( String arg) 
    {
		IJ.register(getClass());
		do
		{
			int res = 0;
			do
			{
				if(!setupDialog())
				{	return;}
				res = validateInput();
				if(res > 1)
				{	IJ.error("bad input! Better luck next time");}
			} while (res > 1);
			
			if(res > 1)
			{	IJ.log("Hex_Symetry: input was autocorrected");}
			
			if(    (current_bondlength != bondlength) 
				|| (current_subframesize != subframesize) )
			{	
				init_transformations();
			}
			apply_transformations();
			
			OutputImg.show();
			OutputImg.updateAndRepaintWindow();
		}
		while(!quit_after_generation);
	
	}
		
	public void apply_transformations()
	{
		if(calculate_median)
		{
			full_histo = new short[subframesize * subframesize][grayvalues]; //reset all
			median = new short[subframesize * subframesize]; //reset all
		}
		else
		{
			full_histo = null;
			median = null;
		}
		int oust = 1;
		ImageStack inputSt = InputImg.getStack();
		ImageStack outputSt = OutputImg.getStack();
		int prog_step = subframeCount/100;
		for(int inst = firstframe; inst <= lastframe; ++inst)
		{		
			short [] inpixels = (short[])inputSt.getPixels(inst);
			for(int mirror = 0; mirror < 2; ++mirror)
			{
				for(int rot = 0; rot < 6; ++rot)
				{
					for(int lp = 0; lp < lpMax; ++lp) // go over all 6-rings in the hex_domain
					{	
						if(!detached && (prog_step != 0) && oust%prog_step == 0)
						{
							IJ.showProgress(oust, subframeCount);
						}
						short[] outpixels = (short[])outputSt.getPixels(oust++);
						for(int hp = 0; hp < hpMax; ++hp) //go over all pixels inside hex_domain
						{
							int ind_i = hexPixels[hp];
							int q_i = ind_i % subframesize;
							int r_i = ind_i / subframesize;
							int ind_o = transform_qr(q_i, r_i, lp, rot, mirror);
							outpixels[ind_o] = inpixels[ind_i];
							if(calculate_median)
							{
								++full_histo[ind_o][inpixels[ind_i]];
							} 
						}
					}
				}
			}
		}
		if(!detached)
		{	IJ.showProgress(1.0);}
		
		
		if(calculate_median)
		{
			final int limit = subframeCount/2;
			for(int hp = 0; hp < hpMax; ++hp) //go over all pixels inside hex_domain
			{
				int ind_i = hexPixels[hp];
				short[] histo = full_histo[ind_i];
				int sum = 0;
				int val = 0;
				while(sum < limit)
				{
					sum += histo[val++];
				}
				median[ind_i] = (short)val;
			}
			if(subtract_median)
			{
				final int area = subframesize * subframesize;
				for(int sl = 1; sl <= subframeCount; ++sl)
				{
					short[] outpixels = (short[])outputSt.getPixels(sl);
					for(int ind = 0; ind < area; ++ind) //25% will be both zero
					{
						outpixels[ind] += (short)(32768 - median[ind]);
					}
				}
			}
			if(show_median)
			{
				ImageStack med_St = ImageStack.create(subframesize,subframesize,1,16);
				med_St.setPixels(median,1);
				ImagePlus med_Img = new ImagePlus(WindowManager.makeUniqueName("MED_" + Outputtitle) , med_St);
				med_Img.show();
			}
			
			
		}
		
		//release the arrays to GarbageCollector
		full_histo = null; 
	}
	
	
	//show Dialog and read all the input
	private boolean setupDialog()
	{
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Hex_Symmetry");
		int[] idArray = WindowManager.getIDList(); // list of all opened images (IDs)
		if (idArray == null)
		{
			IJ.noImage();
			return false;
		}
		
		Outputtitle = WindowManager.makeUniqueName(Outputtitle); // dont overrite previous output
		String[] titleArray = new String[idArray.length]; // titles of opened images
		for (int i = 0; i < idArray.length; ++i)
		{
			titleArray[i] = WindowManager.getImage(idArray[i]).getTitle();
		}
		if (Inputtitle == null || Inputtitle.equals(""))
		{	Inputtitle = titleArray[titleArray.length-1];} //the most recent image
		
		gd.addMessage("The output will be a stack of all possible\n" +
					  "graphene lattice symetry operations.\n" +
					  "The ordering is Tx,Tz,R,M,frames");
		
		gd.addChoice("Source Image/Stack (GRAY16)", titleArray, Inputtitle );
		gd.addNumericField("First frame to process", firstframe, 0);
		gd.addNumericField("Last frame to process (<1 End)", lastframe, 0);
		gd.addNumericField("bondlength",bondlength,0);
		gd.addNumericField("grayvalues (<2 no median)", grayvalues,0);
		gd.addStringField("Output", Outputtitle, 30);
		gd.addCheckbox("calculate median", calculate_median);
		gd.addCheckbox("subtract median", subtract_median);
		gd.addCheckbox("show median", show_median);
		
		gd.addCheckbox("Quit after generation", quit_after_generation);
		gd.showDialog();
		if (gd.wasCanceled())
		{
			return false;
		}
		Inputtitle = gd.getNextChoice();
		firstframe = (int)gd.getNextNumber();
		lastframe = (int)gd.getNextNumber();
		bondlength = (int)gd.getNextNumber();
		grayvalues = (int)gd.getNextNumber();
		Outputtitle = WindowManager.makeUniqueName(gd.getNextString());
		calculate_median = gd.getNextBoolean();
		subtract_median = gd.getNextBoolean();
		show_median = gd.getNextBoolean();
		
		quit_after_generation = gd.getNextBoolean();
		return true;
	}
		
	// 0 .. ok , 1 .. minor issue, 2 .. bad case
	private int validateInput()
	{
		int res = 0; // result code
		
		if(bondlength < 0)
		{
			res = 1;
			bondlength = -bondlength;
		}
		/*
		else if(bondlength == 0)
		{
			IJ.log("bondlength cannot be 0");
			res = 0;
		}
		*/
		if(!detached)
		{
			InputImg = WindowManager.getImage(Inputtitle);
		}
		
		if( (InputImg.getType() != ImagePlus.GRAY16) )
		{
			res = 2;
			IJ.log("Please provide a GRAY16 formated version of " + Inputtitle);
		}
		
		subframesize = InputImg.getWidth(); //actually diameter
		if( (subframesize%2 != 0) ||
			(subframesize != InputImg.getHeight()) )
		{
			IJ.log(Inputtitle + "must be a square with even edge lengths");
			res = 2;
		}
		
		if ( (bondlength != 0) && ( (subframesize/2) % bondlength != 0)	)
		{
			IJ.log("make sure that radius % bondlength == 0");
			res = 2;
		}	 
		hpMax = ( 3 * subframesize * subframesize) / 4; // pixels per hex domain
		if(bondlength != 0)
		{
			lpMax = ( subframesize * subframesize) / ( 4 * bondlength * bondlength ); // number of 6-rings per hex domain
		}
		if(firstframe > lastframe && lastframe > 0)
		{
				int tmp = lastframe;
				lastframe = firstframe;
				firstframe = tmp;
		}

		if(firstframe < 1)
		{
			firstframe = 1;
		}
		
		if( (lastframe < 1) || (lastframe > InputImg.getStackSize()) )
		{
			lastframe = InputImg.getStackSize();
		}
		
		subframeCount = lpMax * 12 * (lastframe - firstframe + 1);
		
		if(grayvalues < 2)
		{
			calculate_median = false;
		}
		
		
		if(!calculate_median && (subtract_median || show_median) )
		{
			subtract_median = false;
			show_median = false;
			IJ.log("cannot subtract or show median from "+ grayvalues +" grayvalues");
			res = 2;
		}
		
		if(calculate_median && !subtract_median && !show_median)
		{
			calculate_median = false;
			res = 1;
		}
		
		
		
		if(res < 2) // no really bad input so far, skip if not needed
		{
			OutputImg = NewImage.createShortImage(Outputtitle, subframesize , subframesize, subframeCount , NewImage.FILL_BLACK);
		}
		return res;
	}
	
	private void init_transformations()
	{
		current_subframesize = subframesize;
		current_bondlength = bondlength;
		
		//prefetch the equivalent lattice sites inside the hex subframe
		final int rH = subframesize/2;
		hexPixels = new int[hpMax];
		latticePoints = new int[lpMax];
		lpX = new int[lpMax];
		lpZ = new int[lpMax];
		
		int hpC = (0);
		int lpC = (0);
		//odd-r -> cube
		Mx = rH - (rH) / 2;
		Mz = rH;
		
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
				int x_data = dx + Mx;
				int z_data = dz + Mz;
				// cube -> odd-r
				int q_data = x_data + (z_data-(z_data&1)) / 2;
				int r_data = z_data;
				int ind = q_data + r_data * subframesize;
				hexPixels[hpC++] = ind;
				
			}
		}
		int dim = (bondlength > 0 )? subframesize/(2*bondlength) : 1;
		for(int stepz = 0 ; stepz < dim; ++stepz)
		{
			for(int stepx = 0 ; stepx < dim; ++stepx)
			{	
				int dx = bondlength*(2*stepx-stepz);
				int dz = bondlength*(2*stepz-stepx);
				int dy = -dx -dz;
				int[] cube = {dx,dy,dz};
				hexagonal_trap(cube,rH);
				dx = cube[0];
				dy = cube[1];
				dz = cube[2];
				
				//translate to center
				int x_data = dx + Mx;
				int z_data = dz + Mz;
				// cube -> odd-r
				int q_data = x_data + (z_data) / 2;
				int ind = q_data + z_data * subframesize; //z_data = r_data
				lpX[lpC] = dx;
				lpZ[lpC] = dz;
				latticePoints[lpC++] = ind; 
			}
		}
	}
	
	private int transform_qr(int q, int r, int lp, int rot, int mirror)
	{
		//odd-r -> cube 
		int x = (q - (r-(r&1))/2 );
		int z = (r);
		
		//offset around center
		int xs = (x - Mx);
		int zs = (z - Mz);
		
		//translation by lattice point
		final int xt = ( lpX[lp] );
		final int zt = ( lpZ[lp] );
		xs += xt;
		zs += zt;
		int ys = (-xs - zs);
		
		//rotation
		if(rot > 0)
		{
			int[] rotxyz = {xs, ys, zs};
			xs = rotxyz[rot % 3];
			ys = rotxyz[(1+rot) % 3];
			if( (rot&1) == 1) //odd number
			{
				xs = -xs;
				ys = -ys;
			}
			zs = (-xs - ys);
		}
		//mirroring
		if( (mirror&1) == 1)
		{
			zs = ys;
			ys = (-xs - zs);
		}

		//periodic hexagon
		int[] cube = {xs, ys, zs};
		hexagonal_trap(cube, Mz);
		
		//undo offset by center
		xs = cube[0] + Mx;
		zs = cube[2] + Mz;
		
		int qt = xs + (zs)/2;
		return (qt + zs * subframesize); //zs = rt
	}
	
	private boolean hexagonal_trap(int[] cube, int N)
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
}
