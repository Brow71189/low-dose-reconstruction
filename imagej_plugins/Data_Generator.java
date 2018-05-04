import ij.*;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.gui.GenericDialog;
import ij.plugin.*;
import ij.io.FileInfo;
import java.util.Random;

public class Data_Generator implements PlugIn {

	int use_molecules = 1; //there can be different independened molecules on every subframe
	int molecules_placed = 0;
	int[] num_molecules = null;
	int subframeCount = 1; //number of subframes to generate
	int bondlength = 4; //C-C distance in hex pixels
	int output_radius = 24; //size of output data
	int input_radius = 0; //size of input data
	int threshold = 10; //1 is good for generated molecules
	int attempts = 10;
	int dim = 1; //periodicity of subframes in hexes
	int num_trans = 1; //number of translations expected in stacks
	int sympersf = 12; //6 rotations and 2 mirror images
	//int mindist = 0; //approx "diameter" of molecules in bondlengths
	String Graphenetitle = WindowManager.getCurrentImage().getTitle(); //the lattice
	String Outtitle = WindowManager.makeUniqueName("Hex_Data.tif"); //the output
	String[] Moleculetitles = null; //the molecules
	double[] density = null; //independent densities of each molecule
	double default_density = 0.05;
	ImagePlus[] MolImgP = null; //the molecul images
	ImagePlus GrapheneImgP; //the lattice images
	ImagePlus OutImgP; //the output image
	Random random = new Random();
	int cx; //hex center of the output
	int cz; //hex center of the output 
	boolean sync_rots_and_mirrors = false; //all molecules and the lattice only differ by translations
	boolean quit_after_generation = false; //quit and dont show the next dialog, maybe good for batchrun
	boolean can_rerun = false;
	boolean lock_dimensions = false;
	boolean slave = false;
	boolean in_place = true;
	
	
	
	
	
	
	public Data_Generator(){} //explicit default constructor
	
	public Data_Generator(String _Outtitle, int _subframeCount, int _bondlength)
	{
		Outtitle = _Outtitle;
		subframeCount = _subframeCount;
		bondlength = _bondlength;
		lock_dimensions = true;
		quit_after_generation = true;
	}
	
	private boolean xorshifter = true; //no noticeable performance cost
    private boolean xorshift_initialized = false;
    private long[] states = new long[16];
	private int p_xor = random.nextInt(16);
	private final double range_H = (double)Long.MAX_VALUE;
	
	private double xorshift1024star() //credits to wikipedia
	{
		double next = 1.0D;
		do
		{
			long s0 = states[ p_xor ];
			long s1 = states[ p_xor = ( (p_xor+1) & 15) ];
			s1 ^= (s1 << 31);
			s1 ^= (s1 >>> 11);
			s0 ^= (s0 >>> 30);											//clear the leading sign bit.
			next = ( (double)( ( ( ( states[p_xor] = s0 ^ s1 ) * 1181783497276652981L ) << 1 ) >>> 1 ) ) / range_H;
			//if(next >= 1.0D || next < 0.0D ) System.out.println("next = " + next);
			
		} 
		while (next >= 1.0); //That should be an extremely rare case, but saftey first
		//The range is [0.0, 1.0[
		return next; 
		
	}
    
    public void init_xorshifter() 
    {
		xorshifter = true;
		//dont reinitialize during batch runs
		if(!xorshift_initialized)
		{
			for(int i = 0; i < states.length; ++i)
			{	//why the hell is that not the default implementation of random.nextLong() !?
				long high = (long)random.nextInt();
				long low = (long)random.nextInt();
				states[i] = ( (high << 32) | low );
			}
			xorshift_initialized = true;
		}
	}
	
	
	
	
	private boolean setupDialog() //show Dialog and read all the input
	{
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog( getClass().getSimpleName() );
		int[] idArray = WindowManager.getIDList(); // list of all opened images (IDs)
		if (idArray == null)
		{
			IJ.noImage();
			return false;
		}
		
		Outtitle = WindowManager.makeUniqueName(Outtitle); // dont overrite previous output
		String[] titleArray = new String[idArray.length]; // titles of opened images
		for (int i = 0; i < idArray.length; ++i)
		{
			titleArray[i] = WindowManager.getImage(idArray[i]).getTitle();
		}
		if (Graphenetitle == null || Graphenetitle.equals(""))
		{	Graphenetitle = titleArray[titleArray.length-1];} //the most recent image
		
		gd.addMessage("Stacks should be arranged by translations > rotations > mirrors\n" +
					  "use Hex_Symmetry to generate them");
		gd.addChoice("Lattice stack (GRAY16)", titleArray, Graphenetitle );
		if(!lock_dimensions) {	gd.addNumericField("bondlength",bondlength,0); }
		//gd.addNumericField("Minimal molecule distance (bondlengths)",mindist,0);
		gd.addNumericField("intensity threshold for collision", threshold, 0);
		gd.addNumericField("attemps for non colldiding placing", attempts, 0);
		for(int i = 0; i < use_molecules; ++i)
		{
			gd.addChoice("Molecule" + i + " stack( GRAY16)", titleArray, Moleculetitles[i] );
			gd.addNumericField("Number of Molecule" + i, num_molecules[i], 0);
			gd.addNumericField("density of Molecule" + i, density[i], 4);
		}
		if(!slave)	{ gd.addStringField("Output", Outtitle, 30);}
		gd.addCheckbox("use xorshiftK* ", xorshifter);
		gd.addNumericField("Number of frames", subframeCount, 0);
		gd.addNumericField("Output_radius", output_radius, 0);
		gd.addCheckbox("replace lattice by molecule", in_place);
		gd.addMessage("This feature totally relies on the correct ordering of the stacks");
		gd.addCheckbox("Synchronize rotations and mirrors", sync_rots_and_mirrors);
		if(molecules_placed != 0)
		{	gd.addMessage("total molecules placed (in last run): " + molecules_placed);}
		gd.addCheckbox("Finish after generation", quit_after_generation);
		gd.showDialog();
		if (gd.wasCanceled())
		{
			return false;
		}
		Graphenetitle = gd.getNextChoice();
		if(!lock_dimensions) {	bondlength = (int)gd.getNextNumber(); }
		//mindist = (int)gd.getNextNumber();
		threshold = (int)gd.getNextNumber();
		attempts = (int)gd.getNextNumber(); 
		for(int i = 0; i < use_molecules; ++i)
		{
			Moleculetitles[i] = gd.getNextChoice();
			num_molecules[i] = (int)gd.getNextNumber();
			density[i] = gd.getNextNumber();
		}
		if(!slave)	{	Outtitle = WindowManager.makeUniqueName(gd.getNextString());}
		xorshifter = gd.getNextBoolean();
		
		subframeCount = (int)gd.getNextNumber();
		output_radius = (int)gd.getNextNumber();
		in_place = gd.getNextBoolean();
		sync_rots_and_mirrors = gd.getNextBoolean();
		quit_after_generation = gd.getNextBoolean();
		//cube hex center for the output
		cx = output_radius - output_radius/2; //this version also works for odd radii
		cz = output_radius;
		return true;
	}
	
	// 0 .. ok , 1 .. minor issue, 2 .. bad case
	private int validateInput()
	{
		int res = 0; // result code
		if (subframeCount < 0)
		{
			res = 1;
			subframeCount = -subframeCount;
		}
		else if(subframeCount == 0)
		{	res = 2;}
		
		if(bondlength < 0)
		{
			res = 1;
			bondlength = -bondlength;
		}
		else if(bondlength == 0)
		{
			IJ.log("bondlength cannot be 0");
			res = 2;
		}
		GrapheneImgP = WindowManager.getImage(Graphenetitle);
		if( (GrapheneImgP.getType() != ImagePlus.GRAY16) )
		{
			res = 2;
			IJ.log("Please provide a GRAY16 formated version of " + Graphenetitle);
		}
		
		input_radius = GrapheneImgP.getWidth(); //actually diameter
		if( (input_radius%(2) != 0) ||
			(input_radius != GrapheneImgP.getHeight()) )
		{
			IJ.log(Graphenetitle + "must be a square with even edge lengths");
			res = 2;
		}
		input_radius = input_radius/2; // diameter -> radius
		if ( (input_radius % bondlength != 0) ||
			 (output_radius % bondlength != 0) ||
			 (output_radius > input_radius)	)
		{
			IJ.log("make sure that radii % bondlength == 0 and/or input_radius >= output_radius");
			res = 2;
		}	 
		dim = (input_radius/bondlength);
		num_trans = dim * dim; //number of hexagons aka eqivalent translations
		sympersf = 12 * num_trans; //total number of symetry elements
		if(GrapheneImgP.getStackSize() != sympersf)
		{
			res = 1;
			sync_rots_and_mirrors = false; //we could only do that with complete stacks
			//IJ.log(Graphenetitle + " appears not to be an exact set of equivalent Images");
		}
		
		for(int i = 0; i < use_molecules; ++i)
		{
			if(density[i] < 0)
			{
				res = 1;
				density[i] = -density[i];
			}
			if(density[i] > 1.0)
			{
				res = 1;
				density[i] = 1.0;
			}
			if( num_molecules[i] < 0)
			{
				res = 1;
				num_molecules[i] = 0;
			}
			
			MolImgP[i] = WindowManager.getImage(Moleculetitles[i]);
			if( (MolImgP[i].getWidth() != 2*input_radius) ||
				(MolImgP[i].getHeight() != 2*input_radius) )
			{
				IJ.log("Dimension missmatch of " + Moleculetitles[i] + " and " + Graphenetitle);
				res = 2;
			}	
			if( (MolImgP[i].getType() != ImagePlus.GRAY16) )
			{
				res = 2;
				IJ.log("Please provide a GRAY16 formated version of " + Moleculetitles[i]);
			}
			if ( MolImgP[i].getStackSize() != sympersf )
			{
				res = 1;
				sync_rots_and_mirrors = false; //we could only do that with complete stacks
				IJ.log(Moleculetitles[i] + " appears not to be an exact set of equivalent Images");
			}
		}
		if(res < 2) // no really bad input so far, skip if not needed
		{	OutImgP = NewImage.createShortImage(Outtitle, 2 * output_radius , 2 * output_radius, subframeCount , NewImage.FILL_BLACK);}
		return res;
	}
	
	
	public void run( String arg) 
    {
		if(arg.equals("Hex_Magic")) 
		{	
			IJ.log("Data_Generator for Hex_Magic");
			slave = true;	
		}
		use_molecules = (int)IJ.getNumber("number of distinct molecules", use_molecules);
		if(use_molecules < 1) 
		{	return;}
		num_molecules = new int[use_molecules];
		for(int i=0; i < use_molecules; ++i)
		{
			num_molecules[i]=1; //default is 1 molecule of each type
		}
		Moleculetitles = new String[use_molecules];
		density = new double[use_molecules];
		for(int i=0; i < use_molecules; ++i) density[i] = default_density;
		MolImgP = new ImagePlus[use_molecules];

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
			{	IJ.log("Data_Generator: input was autocorrected");}
			
			if(quit_after_generation)
			{	//provide some log info if Dialog would not reappear
				IJ.log("lattice: " + Graphenetitle);
				FileInfo graphene_info = GrapheneImgP.getOriginalFileInfo();
				if(graphene_info != null)
				{IJ.log("from file: " + graphene_info.directory + graphene_info.fileName);}
				for(int i = 0; i < use_molecules; ++i)
				{	
					IJ.log("molecule" + i +": num = " + num_molecules[i] + " density = " + density[i] + " " + Moleculetitles[i]);
					FileInfo molecule_info = MolImgP[i].getOriginalFileInfo();
					if(molecule_info != null)
					{IJ.log("from file: " + molecule_info.directory + molecule_info.fileName);}
				}
			}
			generateData();
			if(quit_after_generation)
			{	IJ.log("total molecules placed: " + molecules_placed + " area coverage: " + Math.pow( (double)output_radius/(double)input_radius, 2));}
			
			OutImgP.show();
			OutImgP.updateAndRepaintWindow();
			if(slave)
			{	WindowManager.setTempCurrentImage(OutImgP); }
		}
		while(!quit_after_generation);
		can_rerun = true;
	}
		
	//expects a small "frame" and adds a larger "toadd" to
	//the hexagonal domain of
	public ImagePlus rerun(String data_title, int new_frames)
	{
		if(!can_rerun)
		{
			IJ.log("Data_Generator cannot rerun without initial setup");
			return null;
		}
		subframeCount = new_frames;
		Outtitle = data_title;
		OutImgP = NewImage.createShortImage(Outtitle, 2 * output_radius , 2 * output_radius, subframeCount , NewImage.FILL_BLACK);
		generateData();
		OutImgP.show();
		return OutImgP;
	}
	
	private void addtoframe(short[] frame, short[] toadd)
	{
		int inp_D = 2 * input_radius; //input diameter
		int out_D = 2 * output_radius; //output diameter
		int d_r   = input_radius - output_radius; //assumed to be positive
		//TODO test if explicit hex looping performs better than this
		for( int pxd = 0; pxd < frame.length; ++pxd)
		{
			//normal q,r image coords
			int qo = pxd % out_D;
			int ro = pxd / out_D;
		
			//cube coords with respect to center of output
			int dx = qo - ro/2 - cx;
			int dz = ro - cz;
			int dy = -dx - dz;
			if (dx >= -output_radius && dx <  output_radius &&
					dy >  -output_radius && dy <= output_radius &&
					dz >= -output_radius && dz <  output_radius)
			{
				//we are inside the hex domain of frame
				//apply translation to match the centers of frame and toadd
				int qi = qo + d_r;
				int ri = ro + d_r;
				int pxl = ri * inp_D + qi;
				//unsigned -> signed
				int val = (toadd[pxl] & 0xffff);
				if(in_place)
				{	frame[pxd] = (short) val;}
				else
				{	frame[pxd] += (short) val;}
			}
		}
	}
	
	private int overlapp(short[] frame1, short[] frame2)
	{
		int colliding = 0;
		for( int px = 0; px < frame1.length; ++px)
		{
			if( (frame1[px] > threshold) && (frame2[px] > threshold) )
			{
				++colliding;
			}
		}
		return colliding;
	}



	//choose random frames from the graphene stack
	//pick randomly frames for the molecule stacks
	//the probabilities of each molecule are inependent
	//still the second or third molecule may be rejected
	//if the are to many times to close to earlier molecules
	private void generateData()
	{   
		if(xorshifter) {	init_xorshifter();}
		molecules_placed = 0;
		ImageStack OutSt = OutImgP.getStack();
		int prog_step = subframeCount/100;
		int max_molecules = 0;
		for(int i = 0; i < use_molecules; ++i)
		{	max_molecules += num_molecules[i];}
		
		int[] i_from_mol = new int[max_molecules];	
		{
			int mol = 0;
			for(int i = 0; i < use_molecules; ++i)
			{
				for(int j = 0; j < num_molecules[i]; ++j)
				{
					i_from_mol[mol++] = i;
				}
			}
		}
		int[] perm = new int[max_molecules];
		for(int i = 0; i < max_molecules; ++i) {perm[i] = i;}
		
		for(int sf = 1; sf <= subframeCount; ++sf)
		{
			if( (prog_step > 0) && (sf-1)%prog_step == 0)
			{	IJ.showProgress(sf-1,subframeCount-1);}
			//int grapheneID = 1 + random.nextInt( GrapheneImgP.getStackSize() );
			int maxInt = GrapheneImgP.getStackSize();
			int grapheneID = 1 + ( (int)( (xorshifter ? xorshift1024star() : random.nextDouble()) * (maxInt) ) );
			short[] molecule = null;
			short[] lattice = (short[]) GrapheneImgP.getStack().
				getPixels( grapheneID );
			short[] data = (short[]) OutSt.getPixels(sf); //all black
			
			//TODO use an array to remember all previous loactions and molecules
			boolean occupied = false;
			shuffle(perm);
			for(int sh = 0; sh < max_molecules ; ++sh)
			{
				int mol = perm[sh];
				int i = i_from_mol[mol];
				
				if(((xorshifter ? xorshift1024star() : random.nextDouble()) < density[i]))
				{
					boolean skip = false; 
					int trials = 0;
					do
					{
						//int molT = 1 + random.nextInt(num_trans);
						int molT = 1 + ( (int)( (xorshifter ? xorshift1024star() : random.nextDouble()) * (num_trans) ) );
						int molID = 0;
						if(sync_rots_and_mirrors)
						{	molID = num_trans * ( (grapheneID-1) / num_trans ) + molT;}
						else
						{	
							//molID = 1 + random.nextInt( MolImgP[i].getStackSize() ); 
							int maxmolID = MolImgP[i].getStackSize();
							molID = 1 + ( (int)( (xorshifter ? xorshift1024star() : random.nextDouble()) * (maxmolID) ) );
						}
						molecule = (short[]) MolImgP[i].getStack().getPixels( molID );
						short[] cutmolecule = new short[data.length]; //all zeros
						addtoframe(cutmolecule,molecule);
						if (occupied)
						{
							//test if molecule overlapps with already put molecules
							skip |= (overlapp(cutmolecule, data) > 0);
						}
						if(!skip)
						{
							for(int k=0; k < data.length; ++k)
							{	data[k] += cutmolecule[k];}
							occupied = true;
							++molecules_placed;
						}
					}while(skip && (trials++ < attempts));
				}
			}
			addtoframe(data, lattice);
		}  
		IJ.showProgress(1.0);	
	}

	private void shuffle(int[] perm)
	{
		for(int i = 1; i < perm.length; ++i )
		{
			//int pick = random.nextInt(i+1); // 0..i												
			int pick = (int)( (xorshifter ? xorshift1024star() : random.nextDouble()) * (i+1) );
			/*if(pick < 0)
			{
				IJ.log("Caught a negative pick: " + pick);
			}
			else*/ 
			if( (i != pick) ) 
			{
				int tmp = perm[pick];
				perm[pick] = perm[i];
				perm[i] = tmp;
			} 	
		}
	}

}
