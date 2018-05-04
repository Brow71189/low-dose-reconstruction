import ij.*;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.gui.GenericDialog;
import ij.plugin.*;
import java.util.Random;
import org.apache.commons.math3.distribution.GammaDistribution;

public class Poisson_Noise implements PlugIn {
    ImagePlus imData;
    ImagePlus imNoise;
    double unitsize = 1.0;
    double invunitsize = 1.0;
    double invunitsize2 = 1.0;
    boolean gray8 = false;
    Random random = new Random();

	int firstframe = 1;
	int lastframe = -1;
	int copies = 1;
    boolean filter = false;
    int min_sum = 0;
    int max_sum = -1;
    
    int noiseStackSize = 1;
    int dataStackSize = 1;
    int area = 1;
    int max_level = 2;
    
    GammaDistribution gammaDist = new GammaDistribution( 2.0, 5.0);
    //public boolean skip_values = false;
    
    
    private boolean xorshifter = true; //no noticeable performance cost
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
		for(int i = 0; i < states.length; ++i)
		{	//why the hell is that not the default implementation of random.nextLong() !?
			long high = (long)random.nextInt();
			long low = (long)random.nextInt();
			states[i] = ( (high << 32) | low );
		}
	}
    
    public void run( String arg) 
    {
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Poison_Noise");
		int[] idArray = WindowManager.getIDList(); // list of all opened images (IDs)
		if (idArray == null)
		{
			IJ.noImage();
			return;
		}
		
		String Noisetitle = WindowManager.makeUniqueName("Poison_Noise.tif");
		String[] titleArray = new String[idArray.length]; // titles of opened images
		for (int i = 0; i < idArray.length; ++i)
		{
			titleArray[i] = WindowManager.getImage(idArray[i]).getTitle();
		}
		String Datatitle = WindowManager.getCurrentImage().getTitle();
		if (Datatitle == null || Datatitle.equals(""))
		{	Datatitle = titleArray[titleArray.length-1];} //the most recent image
		gd.addChoice("Source image (GRAY16)", titleArray, Datatitle );
		gd.addNumericField("First frame to process", firstframe, 0);
		gd.addNumericField("Last frame to process (<1 All)", lastframe, 0);
		gd.addNumericField("Number of copies", copies, 0);
		gd.addStringField("Output", Noisetitle, 30);
		//gd.addCheckbox("GRAY8",gray8);
		gd.addNumericField("inv Scaling of mean", unitsize, 3);
		//gd.addCheckbox("draw extra random numbers ", skip_values);
		gd.addCheckbox("use xorshiftK* ", xorshifter);
		gd.addCheckbox("filter by intensity ",filter);
		gd.addNumericField("accepted minimum", min_sum, 0);
		gd.addNumericField("accepted maximum (<0 any)", max_sum, 0);
		
		gd.showDialog();
		if (gd.wasCanceled())
		{
			return;
		}
		Datatitle = gd.getNextChoice();
		firstframe = (int)gd.getNextNumber();
		lastframe = (int)gd.getNextNumber();
		copies = (int)gd.getNextNumber();		
		Noisetitle = WindowManager.makeUniqueName(gd.getNextString());
		//gray8 = gd.getNextBoolean();
		unitsize = gd.getNextNumber();
		invunitsize = 1.0 / unitsize;
		invunitsize2 = 1.0 / ( unitsize * unitsize );
		//skip_values = gd.getNextBoolean();
		xorshifter = gd.getNextBoolean();
		if(xorshifter) {	init_xorshifter();}
		filter = gd.getNextBoolean();
		min_sum = (int)gd.getNextNumber();
		max_sum = (int)gd.getNextNumber();
		imData = WindowManager.getImage(Datatitle);
		if( (imData.getType() != ImagePlus.GRAY16) )
		{
			IJ.error("Please provide a GRAY16 formated version");
			return;
		}
		
		int width = imData.getWidth();
		int height = imData.getHeight();
		dataStackSize = imData.getStackSize();
		area = width*height;
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
		if(lastframe < 1)
		{
			lastframe = dataStackSize;
		}
		noiseStackSize = (1 + lastframe - firstframe ) * copies;
		imNoise = gray8? NewImage.createByteImage(Noisetitle, width , height , noiseStackSize , NewImage.FILL_BLACK)
						:NewImage.createShortImage(Noisetitle, width , height , noiseStackSize , NewImage.FILL_BLACK);
		//for(int i = 0; i < 10; ++i)
		//{	IJ.log("xorshift" + i + "  " + xorshift1024star());}
		
		makeNoise();
		imNoise.show();
    }

	private void makeNoise()
	{
		ImageStack dataStack = imData.getStack();
		ImageStack noiseStack = imNoise.getStack();
		int data_frame = firstframe - 1;
		int noise_frame = 1;
		int progstep = Math.max(noiseStackSize / 100, 1);
		IJ.showProgress(noise_frame, noiseStackSize);
		IJ.showStatus("Poisson_Noise " + noiseStackSize + " Frames");
		short[] noise_pix16 = null;
		byte[] noise_pix8 = null;
		do
		{
			//cycle data between first and lastframe (both inclusive)
			if(++data_frame > lastframe)
			{	data_frame = firstframe;}
			short[] data_pix = (short[])dataStack.getPixels(data_frame);
			if(gray8)
			{	noise_pix8 = (byte[])noiseStack.getPixels(noise_frame);}
			else
			{	noise_pix16 = (short[])noiseStack.getPixels(noise_frame);}
			int sum = 0;
			do
			{
				sum = 0;
				if(gray8)
				{	
					for(int i = 0; i < area; ++i)
					{
						noise_pix8[i] = (byte)poissonValue((double)data_pix[i],0);
						sum += (int)noise_pix8[i];
					}		
				}
				else
				{	
					for(int i = 0; i < area; ++i)
					{
						noise_pix16[i] = (short)poissonValue((double)data_pix[i],0);
						sum += (long)noise_pix16[i];
					}	
				}
			}
			while( filter && ( sum < min_sum || ( (max_sum > -1) && (sum > max_sum) ) ));
			
			if( (++noise_frame % progstep) == 0)
			{
				IJ.showProgress(noise_frame, noiseStackSize);
			}
		} while (noise_frame <= noiseStackSize);
		IJ.showProgress(1.0);
	}
    public void set_unitsize(double new_unitsize)
    {
		invunitsize = 1.0 / new_unitsize;
	}

    /**
     * Algorithm poisson random number (Knuth). While simple, the complexity is
     * linear in λ (the mean).
     *
     * @return a random Poisson-distributed number.
     */
    public int poissonValue(double noiseMean, int level) 
    {
        if(noiseMean == 0)
        {
			return 0;
		}
		else if(noiseMean < 0)
		{
			//out of range
			return -1;
		}
        double limit = Math.exp(-noiseMean*invunitsize);
        if (limit <= 0.0  )
        {
            if(level < max_level)
            {
				//split the big number
				int noiseval = 0;
				int lvl2 = level + 1;
				for(int i=0; i<10; ++i)
				{	noiseval += poissonValue(0.1*noiseMean, lvl2);}
				return noiseval;
			}
			else //out of range
			{	return -1;}
        }
        else
        {
            int k = -1;
            double p = 1.0;
            do 
            {
                ++k;
                p *= (xorshifter ? xorshift1024star() : random.nextDouble());
            } while (p >= limit);
            return k;
        }
    }
    
}
