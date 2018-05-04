import ij.*;
import ij.process.*;
import ij.gui.*;

import java.util.TreeMap;

public class JobData
{
	public double modelmin = 0;
	public double modelmax = 1;
	public double modelmean = 1;
	
	//use these to limit the range of total intensities during applying poison noise
	//-1 disables the respective limit
	public int min_intensity = -1;
	public int max_intensity = -1;
	
	public int datamin = 0;
	public int datamax = 1;
	public long datasum = 1;
	public double datamean = 1.0;
	public double data_bg = 0.0;
	public double data_bg_std = 0.0;
	public long data_num_hexes = 0;
	public double dataStdDev = 1.0;
	public double noise_level = 0.0;
	public double g_theta = 1.0;
	//public int g_k = 1;
	public double g_k = 1.0;
	//private double inv_g_k_fak = 1.0;
	//public double inv_gamma_k = 1.0;
	public double ln_Gk = 0.0;
	
	public double pvaloffset = 0;
	public int use_frames = 0;
	
	public String path = null;
	
	public ImagePlus Idata = null;
	public ImagePlus Ptable = null;
	public long[] dataSumImg = null;
	
	public int threadnum = 0;
	public int[] portions = null;
	public int[] perm = null;
	public long[] histogram = null;
	public int[] zipFrameLen = null;
	public int[] sfmaxval = null;
	public long[] zipImgLen = null;	
	public double[] performances = null;
	
	//Enforce proper use of JobData
	private boolean has_threadnum = false;
	private boolean has_frames = false;
	private boolean has_stats = false;
	private boolean has_portions = false;
	private boolean quick_mode = false;
	
	public void cancel_portions() //needs to be called when perm is changed
	{
		has_portions = false;
	}
	
	public boolean confirm_portions()
	{
		return has_portions;
	}
	
	public void set_threadnum( int workers)
	{
		threadnum = workers;
		portions = new int[workers];
		zipImgLen = new long[workers];
		has_threadnum = true;
		has_stats = false;
	}
	
	public void set_use_frames( int frames)
	{
		use_frames = frames;
		perm = new int[frames];
		for(int i = 0; i<frames; ++i)
		{	perm[i] = i + 1; }
		zipFrameLen = new int[frames];
		sfmaxval = new int[frames];
		has_frames = true;
		has_stats = false;
	}
	
	public boolean fill_portions()
	{
		if(!has_stats || !has_threadnum)
		{
			IJ.error("Cannot fill_portions without stats AND threads");
			return false;
		}
		has_portions = true;
		double dframestart = 0;
		double dframestop = 0;
		int framestop = 0;
		int framestart = 0;
		for (int i = 0; i < threadnum; ++i)
		{
			dframestop = dframestart + ( performances[i] * ((double)use_frames) );
			framestop = (int)(dframestop ); //- 0.0000000001
			if (framestop >= use_frames) { framestop = use_frames-1;}
			if ( (i+1) == threadnum ) { framestop = use_frames-1;}  //last thread gets remaining chunk of data.
			//this thread will receive 1+framestop-framestart frames
			int framecount = 1 + framestop - framestart;
			portions[i] = framecount;
			long zipLen = 0;
			if(!quick_mode)
			{
				for(int frame = framestart; frame <= framestop; ++frame)
				{
					zipLen += zipFrameLen[perm[frame]-1];
				}
			}	
			zipImgLen[i] = zipLen;
			IJ.log("Job" + Hex_Magic.job_id +", rank" + (i) + ": " + framecount + " frames" + (quick_mode?"":"zipped: " + (zipLen/512) + " kB"));
			//set counter for next thread:
			framestart = framestop + 1;
			dframestart = dframestop;
			if(framecount < 1)
			{
				IJ.error("Job" + Hex_Magic.job_id +", rank" + i + " has no data. Quiting Hex_Magic.");
				return false;
			}
		}
		return true;
		//performances = null;
	}
	
	public boolean init_stats(boolean skip_stats)
	{
		if(!has_frames)
		{
			IJ.error("Cannot init stats without data frames");
			return false;
		}
		quick_mode = skip_stats;
		has_stats = true;
		has_portions = false; 
		Poisson_Noise pn = null;
		if(noise_level > 0.0)
		{
			pn = new Poisson_Noise();
			pn.set_unitsize(noise_level);
			pn.init_xorshifter(); //employ custom random xor-shift generator
			if(Hex_Magic.duplicate_Idata)
			{
				ImagePlus IdataClone = Idata.duplicate();
				IdataClone.setTitle(Idata.getTitle());
				IdataClone.show();
			}
		}
		int dmin = Integer.MAX_VALUE;
        int dmax = Integer.MIN_VALUE;
        long dsum = 0;
        long sqrsum = 0;
        long datacount = 0;
        double bg2_sum = 0;
        ImageStack IdataSt = Idata.getStack();
		int modelsize = Hex_Magic.modelsize;//Idata.getWidth(); //must be same as Idata.getHeight()
		int bondlength = Hex_Magic.bondlength;
		int modelarea = modelsize*modelsize;
		data_num_hexes = 0;
		dataSumImg = new long[modelarea];
        int stackSize = use_frames;
        //TreeMap histmap = new TreeMap();
        long histmap[] = skip_stats? null : new long[65536];
        
        for (int sl=1; sl<=stackSize; ++sl)
		{
			int sfmin = Integer.MAX_VALUE;
			int sfmax = Integer.MIN_VALUE;
			int sfnonzero = 0;
			short[] raw_pixels = ( Idata.getType() == ImagePlus.GRAY16 ) ? (short[])IdataSt.getPixels(sl) : null;
			short[] pixels = ( (Idata.getType() != ImagePlus.GRAY16) || (pn != null) ) ? new short[modelarea] : null;
			int noise_sum;
			int watch_dog = 1000;
			ImageProcessor IdataPr = (Idata.getType() == ImagePlus.GRAY16)? null : IdataSt.getProcessor(sl);
			boolean accepted = true;
			//create a short copy of an frame and optionaly apply noisification until intensity is accepted 
			if( (Idata.getType() != ImagePlus.GRAY16) || (pn != null) )
			{
				do //is this loop a a major slow down for pn != null
				{
					if (--watch_dog == 0) {break;}
					
					noise_sum = 0;
					accepted = true;
					for(int i=0; i<pixels.length; ++i)
					{
						int q = i % modelsize;
						int r = i / modelsize;
						if(Hex_Magic.hex_mask[i])//( Hex_Magic.iswithinHex(q,r, modelsize/2) )
						{
							int val = (Idata.getType() == ImagePlus.GRAY16) ? raw_pixels[i] : IdataPr.get(i);
							if( pn != null)	{	val = (int)pn.poissonValue( (double)val, 0);}
							pixels[i] = (short)val;
							noise_sum += val;
						}
					}
					accepted &= ( !( (min_intensity != -1) && ( noise_sum < min_intensity ) ) );
					accepted &= ( !( (max_intensity != -1) && ( noise_sum > max_intensity ) ) );
				} 
				while( !accepted );	
				if(watch_dog < 1)
				{	IJ.log("Warning watch_dog terminated noisifaction of frame " + sl); }	
			}
			else
			{
				pixels = raw_pixels;
			}
			if(pn != null) //update the noisified Idata
			{
				if( Idata.getType() == ImagePlus.GRAY16 )
				{
					IdataSt.setPixels(pixels,sl);
					//Idata.setStack(IdataSt);
				}
				else
				{
					for(int i=0; i<pixels.length; ++i)
					{
						IdataPr.setf( i, (float)pixels[i]);
					}
				}
			}
			
			for (int i=0; i<pixels.length; ++i)
			{
				
				if(Hex_Magic.hex_mask[i])//( Hex_Magic.iswithinHex(q,r, modelsize/2) )
				{
					++datacount;
					
					int val = (int)pixels[i];
					dsum += (long)val;
					if (val>sfmax) sfmax = val;
					if (val<sfmin) sfmin = val;
					if(skip_stats)
					{	continue;}
					++histmap[val];
					/*Long longVal = new Long(val);
					Long longCount = (Long)histmap.get(longVal);
					longCount = new Long( (longCount == null) ? 1 : longCount.longValue() + 1 );
					histmap.put(longVal, longCount);
					*/
					dataSumImg[i] += (long)val;
					
					
					sqrsum += (long)(val*val);
					//if (val>dmax) dmax = val;
					//if (val<dmin) dmin = val;
					
					
					if (val > 0) ++sfnonzero;
					int q = i % modelsize;
					int r = i / modelsize;
					int x = q - r/2 - Hex_Magic.cx;
					int z = r - Hex_Magic.cz;
					int y = -x -z;
					int[] cube = {x,y,z};
					Hex_Magic.hexagonaltrap(cube, bondlength);
					if( (cube[0] == 0) && (cube[1] == 0) ) //no need to also test cube[2]
					{ //This pixel is at the dark center of a graphene hexagon
						++data_num_hexes;
						data_bg += (double)val;
						bg2_sum += Math.pow((double)val,2);
					}
				}
			}
			
			sfmaxval[sl-1] = sfmax;
			zipFrameLen[sl-1]= 2*(1+sfmax-sfmin) + sfnonzero; //  2* (1+sfmax-sfmin)+sfnonzero
			if (sfmax>dmax) dmax = sfmax;
			if (sfmin<dmin) dmin = sfmin;
		}
		datasum = dsum;
        datamin = dmin;
        datamax = dmax;
        if( Idata.getType() == ImagePlus.GRAY16 )
		{
			Idata.setStack(IdataSt);
		}
        if(skip_stats)
        {	
			return true;
		}
		histogram = new long[1 + dmax - dmin];
		for(int val = dmin; val <= dmax; ++val)
		{
			/*Long longVal = new Long(val);
			Long longCount = (Long)histmap.get(longVal);
			long count = (longCount == null) ? 0 : longCount.intValue();*/
			long count = histmap[val];
			histogram[val-dmin] = count;
			//System.out.println("Job" + job + "\tval: " + val + "\tcts.: " + count);
		}
		data_bg /= (double)data_num_hexes;
		double data_bg_var = bg2_sum/(double)data_num_hexes - Math.pow(data_bg,2);
		data_bg_std = Math.sqrt(data_bg_var);
		g_k = 1.0;
	
		g_theta = data_bg;
		//inv_gamma_k = 1.0;
        if(data_bg > 0.0)
        {
			g_theta = data_bg_var/data_bg;
			g_k = Math.pow(data_bg,2)/data_bg_var;
			ln_Gk = lnGamma(g_k);
		}
        
       
        
        
        double mean = (double)dsum/(double)datacount;
        datamean = mean;
        dataStdDev = Math.pow((double)sqrsum/(double)datacount - mean*mean, 0.5);
/*        
        if(pn != null)
        {
			Idata.updateAndRepaintWindow();
			Idata.show();
		}
*/		
        return true;
	}
	
	//////////
	// can we define a reasonable offset for the prob. table?
	// mean value from data should be equal to mean value from model!
	// init model with mean value (+ small noise)
	// and adjust ptable so that P(Data|mean value)=1.
	public void autopvaloffset(double new_pvaloffset, boolean auto_offset)
	{
		pvaloffset = new_pvaloffset; //will be overwritten by auto_offset
		if(auto_offset)
		{
			double apvaloffset = 0.0;
			double inv_num = 1/( (double)((long)use_frames * (long)Idata.getWidth() * (long)Idata.getHeight() * 0.75) );//hexes cover 3/4 of models
			for (int val = datamax; val >= datamin; --val) //logprobs for biggest datavalues are typically the smallest
			{	apvaloffset += (double)histogram[val-datamin] * poissonlogprob(val,datamean);} //These two should be like x and ln(x)
			pvaloffset = apvaloffset * inv_num;
		}
	}
	
	public void createPtable(int modelGraylevels)
	{
		Ptable = NewImage.createFloatImage("Prob. lookup table", 1+datamax-datamin  , modelGraylevels , 1, NewImage.FILL_BLACK);
		ImageProcessor Ptablep = Ptable.getProcessor();
		double pvalfloor= - 60.0;
		double fchecksum = 0.0;
		for (int datav=datamin; datav<(1+datamax-datamin); ++datav) //calculate Ptable
		{
			for (int modelv=0; modelv<modelGraylevels; ++modelv)
			{
				double mval = modelmin+(modelmax-modelmin)*(((double)modelv)/((double)modelGraylevels));
				double pval = 0.0;
				if ((mval == 0) && (datav > 0)  ) //&& (data_bg == 0.0)
				{   // totally impossible with experimental data 
					//hack probabilities for modelv == 0 are a power of modelv == 1
					double mval2 = modelmax / ( (double)modelGraylevels );//modelmin must be zero in that case
					pval = 1.2*poissonlogprob(datav,mval2);//(jointlogprob(datamin + datav, mval2, data_bg));
				}
				else
				{	pval = poissonlogprob(datav,mval);}//(jointlogprob(datamin + datav, mval, data_bg));}
				pval -= pvaloffset;
				if (pval<pvalfloor) //a regularization.
				{	pval=pvalfloor*(1.0d+Math.tanh((pval-pvalfloor)/pvalfloor));}
				Ptablep.putPixelValue(datav,modelv,(float)pval);
				fchecksum += pval;
			}
		}
		IJ.log("Job" + (Hex_Magic.job_id + 1) + " Ptable, offset: " + Hex_Magic.format.format(pvaloffset) +
		", floor: " + Hex_Magic.format.format(pvalfloor) + ", checksum: " + Hex_Magic.format.format(fchecksum) );		
		IJ.log("Bg Gamma dist: theta: " + g_theta + "k: " + g_k );
		Ptable.show();
	}

	
	private double jointlogprob(int k, double r, double bg)
    {
        double jointProb = 0.0;
        for(int k1 = 0; k1 <= k; ++k1)
        {	//k1 counts on background
		  	//k2 counts on model
			int k2 = k - k1;
			//Change this line for other distribution functions of bg and model
			//a little useless for two times the same divisible PDF
			jointProb += poissonprob(k1, bg) * poissonprob(k2, r);
		}
		return Math.log(jointProb);
    }
	
	private double poissonprob(int k, double r)
	{
		long fakk = 1;
		int kk = k;
		while (kk > 1) {	fakk *= kk--;}
		return Math.pow(r,k)*Math.exp(-r)/(double)fakk;
		//return result;	
	}
	
	private double poissonlogprob(int x, double avg)
	{
		if( avg <= 0.0)
		{	return (x==0)?1.0:0.0;}
		if( x == 0)
		{	return Math.exp(-avg);}
		int kk = x;
		double f0 = Math.log(avg);
		double f = (double)x*f0 - avg;
		while (kk > 1)	{	f -= Math.log( (double)(kk--) );}
		return f;
	}
	
	
	
	
	private double expprob(double k, double r)
	{
		return Math.exp(-k/r)/r;
	}
	
	
	private double gaussprob(int k, double r)
	{
		//Normalization is not exact but not need at all anyways
		return 	Math.exp( -0.5*Math.pow( ((double)k-r)/data_bg_std , 2 )  ) * 
				Math.pow( 2*Math.PI, -0.5 ) / data_bg_std;
	}
	
	
	private double pointlogprob(int k, double r)
    {
        if ( (k==0) ) return -r;
        return k*Math.log(r)-logGamma2(k)-r;
    }

    private double logGamma2(int k)
    {
        double lnkfak = Math.log((double)k);
        while(--k > 0)
        {	lnkfak += Math.log((double)k);}
        return lnkfak;
    }
    
	private double lnGamma(double x) 
	{
		// Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
		// Communications of the Association for Computing Machinery, 9:684
		double f = 0.0, z;
		if (x < 7) 
		{
			f = 1;
			z = x - 1;
			while (++z < 7)
			{
				f *= z;
			}
			x = z;
			f = -Math.log(f);
		}
		z = 1 / (x * x);
		return	f + (x - 0.5) * Math.log(x) - x + 0.918938533204673 +
				( ( ( -0.000595238095238  * z + 0.000793650793651 ) * z 
					  -0.002777777777778) * z + 0.083333333333333 ) / x;
	}      
}
