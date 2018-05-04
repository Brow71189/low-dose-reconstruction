import ij.*;
import ij.IJ;
import ij.gui.*;
import ij.process.*;
import ij.text.TextWindow;
import ij.text.TextPanel;
import ij.measure.ResultsTable;
import ij.plugin.filter.GaussianBlur;
import java.util.Arrays;
 
public class HexAnalyzer
{
	boolean interactive_mode = false;
	boolean create_wedge = false;
	boolean create_PDF_maps = false;
	boolean create_prob_map = false;
	boolean create_table = false;
	boolean silent_mode = false;
	boolean interpolate = false;
	boolean virtual_smoothing = false;
	boolean plot_std_vs_avg = false;
	boolean keep_step = true;
	boolean sym_Histo = true;
	String fitname = "linear";
	int impDepth = 1;
	int modelsize = 64;
	int modelarea = modelsize * modelsize;
	int modelnum = 0;
	int bondlength = 4;
	//FIXME hardcoded hack for bondlength = 4
	int[] perm4 = {0,1,3,2,6,4,8,7,5};
	int max_val = 0;
	int min_val = 0;
	int real_min_val = 0;
	int max_pos = 1;
	int min_pos = 1;
	double model_min = 0;
	double model_max = 0;
	double data_mean = 0;
	double pvaloffset = 0.0;
	double pvalfloor = - 3.0;
	int minfactor = 3;
	
	int histwidth = max_val-min_val+1;
	public int min_histwidth = 0; //auto
	int histheight = 1;
	String dataTitle = null;
	ImagePlus Idata = null;
	ImageStack IdataSt = null;
	String histTitle = "<create>";
	ImagePlus histImg = null;
	ImageStack histImgSt = null;
	String[] probArray = {"<none>","PoissonF","GammaF","GaussF","Exp","Histo","invGammaF","invBetaF","Rubber","LogNormF"}; //TODO add "All"
	int probChoice = 8;
	String probTitle = "Rubber";
	ImagePlus probImg = null;
	ImageStack probImgSt = null;
	String probMapTitle = "Prob_Map.tif";
	ResultsTable table = null;
	ResultsTable table2 = null;
	
	int output_level = 1;
	int modelGraylevels = 128; //reasonable default
	double modelBrightness = -4; //good default for defects in graphene
	double brightness_offset = 3; //barely sufficient for void in 8 ring
	double gb_sigma = 1.5; //minimal smoothing
	
	int[] orderedYs = null;
	long[] linesums = null;
	double[] rel_Ws = null; 
	double[] areas = null;
	double[] avgs = null;
	double[] modes = null;
	double[] peak_vals = null;
	double[] lnAvgs = null;
	double[] vars = null;
	double[] stds = null;
	double[] skews = null;
	double[] kurts = null;
	double[] alphas = null;
	double[] betas = null;
	double[] iGNs = null;
	
	//inverse Beta Distribution
	double[] iBalphas = null;
	double[] iBbetas = null;
	double[] iBNs = null;
	
	//Gamma Distribution k,Theta x0
	double[] thetas = null;
	double[] lnThetas = null;
	double[] ks = null;
	double[] lnGks = null;
	
	//log normal distribution
	double[] LNmus = null;
	double[] LNsigs = null;
	double[] LNvars = null;
	
	double[][] matches = null;
	
	double[] x0Gs = null;
	double[] x0iGs = null;
	double[] x0Ps = null;
	double[] x0iBs = null;
	double[] x0LNs = null;
	
	boolean[] hex_mask = null;//that should speed up things a lot
	
	double x0PM = 0.0;
	double x0GM = 0.0;
	double x0iGM = 0.0;
	double x0iBM = 0.0;
	double x0LN = 0.0;
	
	short[] unique = null;
	
	int rH = modelsize/2;
	int cq = modelsize/2;
	int cr = modelsize/2;
	int cx = cq - cr / 2;
	int cz = cr;
	int cy = -cx - cz;
	
	int flags = 0; //0	NOTHING
	
	String called_by = "";
	
	public OutlierCore outlierCore = null;
	private CurveFitter stdFitter = new CurveFitter();
	
	HexAnalyzer(){}; //Illegal Access Exception
	
	private void setup_histImg()
	{
		if(histTitle.equals("<create>"))
		{
			IJ.showStatus("scanning " + Idata.getTitle());
			IJ.showProgress(0);
			boolean first = true;
			int progstep = impDepth/100;
			if(progstep == 0)
			{	progstep = 1;}
			
			for(int sl = 1; sl <= impDepth; ++sl)
			{
				short[] model = (short[]) IdataSt.getPixels(sl);
				for(int i=0; i < modelarea; ++i)
				{
					if( hex_mask[i] )
					{
						int val = (int)(model[i] & 0xffff);
						if(!first)
						{	
							if( val < min_val)
							{	min_val = val; min_pos = sl; }
							if( val > max_val)
							{	max_val = val; max_pos = sl; }
						}
						else 
						{
							first = false;
							min_val = val;min_pos = sl;
							max_val = val;max_pos = sl;
							//IJ.log("first value: " + val + "@" + q + "," + r);
						}	
					}
				}
				if(sl%progstep == 0)
				{	
					//IJ.showStatus("scanning " + Idata.getTitle());
					IJ.showProgress(sl,impDepth+1);
				}
			}
			real_min_val = min_val;
			histwidth = max_val-min_val+1;
			if(min_histwidth > 0)
			{	
				if(max_val >= min_histwidth)
				{
					IJ.log("HexAnalyzer: WARNING max value " + max_val + " exceeds fixed histogramm width " + histwidth);
				}
				else
				{
					histwidth = min_histwidth;
					max_val = histwidth - 1;
					min_val = 0;
				}
			}
			
			
			if(sym_Histo)
			{
				histheight = 0; //number of places on irreducible wedge
				for(int m = 0; m <= bondlength; ++m)
				{	histheight += (1+m/2);} 
			}
			else
			{
				histheight = 3*bondlength*bondlength;
			}
			histImg = NewImage.createFloatImage(WindowManager.makeUniqueName("Hex_Histogramm.tif"),histwidth, histheight, 1, NewImage.FILL_BLACK);
			histImgSt = histImg.getStack();
			if(!silent_mode) {	histImg.show();}
			IJ.showProgress(1.0);
		}
		else
		{
			histheight = histImg.getHeight();
			histwidth = histImg.getWidth();
			max_val = min_val + histwidth - 1;
		}
		
		orderedYs = new int[histheight];
		linesums = new long[histheight];
		rel_Ws = new double[histheight];
		areas = new double[histheight];
		avgs = new double[histheight];
		modes = new double[histheight];
		peak_vals = new double[histheight];
		lnAvgs = new double[histheight];
		vars = new double[histheight];
		
		stds = new double[histheight];
		
		skews = new double[histheight];
		kurts = new double[histheight];
		alphas = new double[histheight];
		
		betas = new double[histheight];
		
		iGNs = new double[histheight];
		
		iBalphas = new double[histheight];
		iBbetas = new double[histheight];
		iBNs = new double[histheight];
		thetas = new double[histheight];
		lnThetas = new double[histheight];
		LNmus = new double[histheight];
		LNsigs = new double[histheight];
		LNvars = new double[histheight];
		x0Gs = new double[histheight];
		x0iGs = new double[histheight];
		x0Ps = new double[histheight];
		x0iBs = new double[histheight];
		x0LNs = new double[histheight];
		ks = new double[histheight];
		lnGks = new double[histheight];
		matches  = new double[histheight][14];
	
		int total_area = 0;
		for(int i=0; i < modelarea; ++i)
		{
			if( hex_mask[i] )
			{
				++total_area;
				int y = unique[i];
				rel_Ws[y] += 1.0;
			}
		}
		for(int y = 0; y < histheight; ++y)
		{	
			rel_Ws[y]/=(double)total_area;	
			orderedYs[y] = y;
		}
	}
	
	private int get_unique_pixel(int q, int r)
	{
		int px = q - r / 2;
		int pz = r;
		int py = -px - pz;
		//offset hex in cube coords
		int dx = px - cx;
		int dy = py - cy;
		int dz = pz - cz;
		int[] cube = {dx, dy, dz};
		
		hexagonaltrap( cube, bondlength);
		int ax = Math.abs(cube[0]);
		int ay = Math.abs(cube[1]);
		int az = Math.abs(cube[2]);
		if(sym_Histo)
		{
			//biggest and smallest absolut value determine position on irreducible wedge
			int mm = ax;
			int nn = ax;
			if(ay > mm) {	mm = ay;}
			if(az > mm) {	mm = az;} 
			if(ay < nn) {	nn = ay;}
			if(az < nn) {	nn = az;}
			switch (mm)
			{
				case 0:
					return 0;
				case 1:
					return 1; 
				case 2:
					return 2 + nn; //nn .. 0,1
				case 3:
					return 4 + nn; //nn .. 0,1
				case 4:
					return 6 + nn; //nn .. 0,1,2
				default: //slower version for bondlength > 4
				int col = 6;
				int m = 4;
				while( m < mm)
				{
					int dn = 1 + m/2;
					col += dn;
					++m;
				}
				return col + nn;
			}
		}
		else
		{
			final int x = cube[0];
			final int y = cube[1];
			final int z = cube[2];
			
			final int d = (ax+ay+az)/2;
			int interior= 0;
			if(d>0)
			{
				interior = 1;
				for(int dd = 1;dd < d; ++dd)
				{
					interior += 6*dd;
				} 
			
				if( x == d )
				{
					interior -= z; 
				}
				else if(z == -d)
				{
					interior += ( d<bondlength ? d : -1 );
					interior += y;
				}
				else if(y == d)
				{
					interior += ( d<bondlength ? 2*d: bondlength-1);
					interior -= x;
				}
				else if(x == -d)
				{
					interior +=  (d<bondlength ? 3*d: 2*bondlength-1);
					interior += z;
				}
				else if(z == d)
				{
					interior += (4*d);
					interior -= y;
				}
				else if(y == -d)
				{
					interior += (5*d);
					interior += x;
				}
			}
			return interior%(3*bondlength*bondlength);
		}
	} 
	
	private void populate_histImg()
	{
		IJ.showStatus("populating " + histImg.getTitle());
		IJ.showProgress(0);
		if(silent_mode) {	virtual_smoothing = false;}
		
		if(virtual_smoothing)
		{
			if(outlierCore == null)
			{	outlierCore = new OutlierCore(flags);}
			outlierCore.create_report = true;
			virtual_smoothing = (outlierCore.showDialog(Idata,"virtual Outlier Masking") == flags);
		}
		
		float[] histpix = (float[])histImgSt.getPixels(1);
		long[] histraw = new long[ histpix.length ]; 
		long total_pixels=0;
		//if (histpix.length != histwidth * histheight) IJ.log("Dimensions of histImg pixels dont match");
		int progstep = impDepth/100;
		if(progstep == 0)
		{	progstep = 1;}
		for(int sl = 1; sl <= impDepth; ++sl)
		{
			short[] pixels = (short[]) IdataSt.getPixels(sl);
			if(virtual_smoothing)
			{
				pixels = Arrays.copyOf(pixels,pixels.length);
				outlierCore.run(pixels);
			}
			for(int ind=0; ind < modelarea; ++ind)
			{
				if( hex_mask[ind] )
				{
					++total_pixels;
					int row = (int)(pixels[ind] & 0xffff) - min_val;
					//if( row >= histwidth) IJ.log( "row exceeds histwidth " + row);
					int col = unique[ind];//perm4[get_unique_pixel(q, r)];
					//if( col >= histheight ) IJ.log( "col exceeds histheight " + col );
					int hist_ind = (row + col * histwidth);
					//if (hist_ind >= histpix.length) IJ.log ("hist_ind: " + hist_ind + " exceeds " + histpix.length);
					//if (hist_ind < 0 ) IJ.log ("hist_ind below 0: " + hist_ind);
					//IJ.log("hist_ind: " + hist_ind + "  row: " + row + "  col: " + col);
					++histraw[hist_ind]; 
				}
			}
			if(sl%progstep == 0)
			{	
				//IJ.showStatus("populating " + histImg.getTitle());
				IJ.showProgress(sl,impDepth+1);
			}
		}
		double total_area = (double) total_pixels;
		for(int y = 0; y < histheight; ++y) 
		{	
			long linesum = 0L;
			for(int x = 0; x < histwidth; ++x) 
			{	linesum += histraw[x + y * histwidth];	}
			
			linesums[y] = linesum;
			if(linesum > 0)
			{
				double area = (double)linesum;
				for(int x = 0; x < histwidth; ++x)
				{
					histpix[x + y * histwidth] = (float)
						( ( (double)histraw[x + y * histwidth]) / area);	
				}
			} 
		}
		if(!silent_mode) {	histImg.updateAndRepaintWindow();}
		IJ.showProgress(1.0);
	}
	
	private void show_inequivalent_pixels()
	{
		//useless in slient mode or without unique pixels
		if(silent_mode || unique == null)
		{	return;}
		ImagePlus ineq = NewImage.createShortImage(WindowManager.makeUniqueName("Equivalent_sites.tif"),modelsize, modelsize, 1, NewImage.FILL_BLACK);
		ImageStack ineqSt = ineq.getStack();
		ineqSt.setPixels(unique,1);
		ineq.setStack(ineqSt);
		ineq.show();
	}
	
	private void make_unique_pixels()
	{
		unique = new short[modelarea];
		for(int i=0; i < modelarea; ++i)
		{
			if( hex_mask[i] )
			//values with +1 match rownumbers in ResultsTable
			{	
				int q = i % modelsize;
				int r = i / modelsize;
				int id = get_unique_pixel(q, r);
				if( sym_Histo && (bondlength == 4) )
				{	id = perm4[id];}
				unique[i] = (short) id;		
			} 
		}
	}
	
	public ImagePlus make_models(String mtitle, int mnum)
	{
		ImagePlus model = NewImage.createShortImage(WindowManager.makeUniqueName(mtitle),modelsize, modelsize, mnum, NewImage.FILL_BLACK);
		ImageStack modelSt = model.getStack();
		short[] pixels = (short[])modelSt.getPixels(1);
		for(int i=0; i < modelarea; ++i)
		{
			if( hex_mask[i] )
			{	
				pixels[i] = (short)((modelGraylevels*(avgs[unique[i]] - model_min))/(model_max-model_min)  + 0.499);
			}
			else
			{	pixels[i] = 0;}
		}
		for(int m = 2; m <= modelnum; ++m)
		{
			short[] pixm = (short[])modelSt.getPixels(m);
			for(int i=0; i < modelarea; ++i)
			{	pixm[i] = pixels[i];}
		}
		if(!silent_mode) {	model.show();} 
		return model;
	}
	
	public ImagePlus make_asym_models(String mtitle, int mnum)
	{
		ImagePlus model = NewImage.createShortImage(WindowManager.makeUniqueName(mtitle),modelsize, modelsize, mnum, NewImage.FILL_BLACK);
		ImageStack modelSt = model.getStack();
		short[] mpixels = (short[])modelSt.getPixels(1);
		long[] sumpix = new long[modelsize*modelsize];
		
		for(int sl = 1; sl <= impDepth; ++sl)
		{
			short[] pixels = (short[]) IdataSt.getPixels(sl);
			if(virtual_smoothing)
			{
				pixels = Arrays.copyOf(pixels,pixels.length);
				outlierCore.run(pixels);
			}
			for(int i = 0; i < pixels.length; ++i)
			{
				sumpix[i] += pixels[i];
			}
		}
		
		for(int i=0; i < modelarea; ++i)
		{
			if( hex_mask[i] )
			{	
				double avg = ((double)sumpix[i])/impDepth;
				double mval =(modelGraylevels-1)*(avg-model_min)/(model_max-model_min);
				mpixels[i] = (short)(mval+0.5);	
			}
			else
			{	mpixels[i] = 0;}
		}

		for(int m = 2; m <= modelnum; ++m)
		{
			short[] pixm = (short[])modelSt.getPixels(m);
			for(int i=0; i < modelarea; ++i)
			{	pixm[i] = mpixels[i];}
		}
		if(!silent_mode) {	model.show();} 
		return model;
	}
	
	
	
	private void do_statistics()
	{
		data_mean = 0.0;
		float[] histpix = (float[])histImgSt.getPixels(1);
		for(int y = 0; y < histheight; ++y) 
		{
			double sum = 0.0;
			double sum2 = 0.0;
			double sum3 = 0.0;
			double sum4 = 0.0;
			double weight = 0.0;
			double mode = 0.0;
			double max = 0.0;
			boolean first = true;
			//first pass for mean
			for(int x = 0; x < histwidth; ++x)
			{
				int xp = x + min_val;// + data_offset;
				double w = (double)histpix[x + y * histwidth];
				if(first || (w > max) )
				{
					first = false;
					max = w;
					mode = xp;
				}
				weight += w;
				sum += (double)xp * w;
				
			}
			double avg = sum/weight;
			double avg2 = Math.pow(avg,2);
			//second pass for variance and standard deviation
			for(int x = 0; x < histwidth; ++x)
			{
				int xp = x + min_val;// + data_offset;
				double w = (double)histpix[x + y * histwidth];
				sum2 += Math.pow(((double)xp-avg),2) * w;
			}
			double variance = (sum2/weight);
			double std = Math.sqrt(variance);
			//third pass for skewness
			for(int x = 0; x < histwidth; ++x)
			{
				int xp = x + min_val;// + data_offset;
				double w = (double)histpix[x + y * histwidth];
				sum3 += Math.pow(((double)xp-avg)/std,3) * w;
				sum4 += Math.pow(((double)xp-avg)/std,4) * w;
				
			}
			double skew = sum3/weight;
			double kurtosis = sum4/weight;
			//Gamma Distribution
			double k = Math.pow(2.0/skew,2);
			double theta = std*skew/2.0;
			double x0G = k*theta - avg;
			double ln_Gk = lnGamma(k);
			//Inverse Gamma Distribution
			double skew2 = Math.pow(skew,2);
			//a=skew2
			double b = (6*skew2 + 4);
			double b2 = Math.pow(b,2);
			double c = 9*skew2-8;
			double alpha = ( b + Math.sqrt( b2 - 4*skew2*c) ) / 
							(2*skew2);
			double beta = Math.sqrt(alpha-2)*(alpha-1)*std;				
			double iGN =  Math.exp(alpha * Math.log(beta) - lnGamma(alpha));
			double x0iG = Math.sqrt(alpha-2)*std - avg;
			//Poisson Distribution
			double x0P = variance-avg;
			//inverse Beta Distribution
			double sk2var = Math.pow(skew,2)*variance;
			double skstd = skew * std;
			//a=sk2var
			b = 6*sk2var - 8*skstd - 16*variance;
			b2 = Math.pow(b,2);
			c = 12 + 32*variance + 24 * skstd;
			double iBbeta = ( b + Math.sqrt( b2 + 4*sk2var*c) ) / 
							(2*sk2var);
			double iBalpha = 0.5*(iBbeta-1)*( 1 + Math.sqrt(1+4*(iBbeta-2)*variance) );
			double iBN = lnEulerBeta(iBalpha,iBbeta);
			double x0iB = iBalpha/(iBbeta-1) - avg;
			//lognormal Distribution
			b = Math.pow(0.5*(skew*Math.pow(4 + skew2,0.5) + skew2 + 2),1.0/3.0);
			double u = (b*b-b+1)/b;
			double LNvar = Math.log(u);
			double LNsig = Math.pow(LNvar,0.5);
			double LNmu = 0.5*Math.log(variance/(Math.pow(u,2) - u));
			double x0LN = Math.exp(LNmu+0.5*LNvar) - avg;
			//Outdated version for Gamma Distribution
			//double theta = variance/avg;
			//double k = avg2/variance;
			
			// Descriptive statistics
			areas[y] = weight;
			avgs[y] = avg;
			modes[y] = mode;
			peak_vals[y] = max;
			lnAvgs[y] = Math.log(avg);
			vars[y] = variance;
			stds[y] = std;
			skews[y] = skew;
			kurts[y] = kurtosis;
			//Gamma Distribution
			thetas[y] = theta; //scale
			lnThetas[y] = Math.log( theta );
			ks[y] = k; //shape
			
			lnGks[y] = ln_Gk;
			x0Gs[y] = x0G;
			//Inverse Gamma Distributin
			alphas[y] = alpha;
			betas[y] = beta;
			iGNs[y] = iGN;
			x0iGs[y] = x0iG; 
			//Poisson Distribution
			x0Ps[y] = x0P;
			//Inverse Beta Distribution
			iBalphas[y] = iBalpha;
			iBbetas[y] = iBbeta;
			iBNs[y] = iBN;
			x0iBs[y] = x0iB;
			//lognormal Distribution
			LNsigs[y] = LNsig;
			LNvars[y] = LNvar;
			LNmus[y] = LNmu;
			x0LNs[y] = x0LN;
			
			data_mean += avg * rel_Ws[y];
			
		}
		
		x0PM  = 0.0;
		x0GM  = 0.0;
		x0PM  = 0.0;
		x0iGM = 0.0;
		x0iBM = 0.0;
		x0LN  = 0.0;
		
		for(int y = 0; y < histheight; ++y) 
		{
			double w = 1.0/histheight;//rel_Ws[y];
			x0PM += w * x0Ps[y];
			x0GM += w * x0Gs[y];
			x0iGM += w * x0iGs[y];
			x0iBM += w * x0iBs[y];
			x0LN += w * x0LNs[y];
		}
		stdFitter.set_pts(avgs,stds);
		if(!fitname.equals("linear") && !fitname.equals("auto"))
		{	
			stdFitter.fit("linear");
			System.out.println(stdFitter.get_formula());
		}
		stdFitter.fit(fitname);
		
		if(fitname.equals("auto"))
		{
			IJ.log("best std vs. avg: " + stdFitter.get_model() + " rms: " + stdFitter.get_merit() );
			int last_mode = stdFitter.get_all_merits().length;
			for(int m = 0; m < last_mode; ++m)
			{
				System.out.println(stdFitter.get_formula(m));
			}
		}
		else
		{
			System.out.println(stdFitter.get_formula());
		}
		//lin_reg();
		//IJ.log("intercept: " + intercept1 + "  slope: " + slope1);
		order_by_key(avgs,orderedYs);
		//for(int y0 = 0; y0 < histheight; ++y0)
		//{	IJ.log("perm4[" + y0 + "] =  " + orderedYs[y0]);}
		
		if(create_table)
		{
			table = new ResultsTable();
			for(int y0 = 0; y0 < histheight; ++y0)
			{	
				int y = orderedYs[y0];
				//table.setValue("mode", y0, mode);
				table.setValue("row", y0, y);
				table.setValue("avg.", y0, avgs[y]);
				table.setValue("std.", y0, stds[y]);
				table.setValue("var.", y0, vars[y]);
				table.setValue("skew", y0, skews[y]);
				table.setValue("kurt.", y0, kurts[y]);
				table.setValue("x0P", y0, x0Ps[y]);
				table.setValue("x0G", y0, x0Gs[y]);
				table.setValue("k", y0, ks[y]);
				table.setValue("theta", y0, thetas[y]);
				table.setValue("x0iG", y0, x0iGs[y]);
				table.setValue("iGalpha", y0, alphas[y]);
				table.setValue("iGbeta", y0, betas[y]);
				table.setValue("x0iB", y0, x0iBs[y]);
				table.setValue("iBalpha", y0, iBalphas[y]);
				table.setValue("iBbeta", y0, iBbetas[y]);
				table.setValue("x0LN", y0, x0LNs[y]);
				table.setValue("LNmu", y0, LNmus[y]);
				table.setValue("LNsig", y0, LNsigs[y]);	
				//table.setValue("lnEuB", y0, iBN);			
			}
			table.showRowNumbers(true);
			table.show("Stats " + (Idata!=null?Idata.getTitle():histImg.getTitle()));
		}
		
		if(plot_std_vs_avg)
		{
			double[] xreg = new double[histwidth];
			double[] yreg = new double[histwidth];
			for(int i=0; i < histwidth; ++i)
			{
				xreg[i] = min_val + i;
				yreg[i] = stdFitter.get_val(xreg[i]); 	
			}
			String s = stdFitter.get_model();
			double[] ylin = new double[histwidth];
			for(int i=0; i < histwidth; ++i)
			{
				ylin[i] = stdFitter.get_val(1,xreg[i]); 	//1 .. "linear"
			}
			
			Plot plot = new Plot( dataTitle+"_"+s , "avg", "std");
			
			plot.setColor("BLACK"); //linear reference sets display range
			plot.addPoints(xreg, ylin, Plot.DEFAULT_FLAGS & Plot.LINE);
			//plot.show();
			plot.setColor("RED"); //red fit function on top
			plot.addPoints(xreg, yreg, Plot.DEFAULT_FLAGS & Plot.LINE);
			//plot.show();
			plot.setColor("BLUE"); //lastly the data points
			plot.addPoints(avgs,stds, Plot.DEFAULT_FLAGS & Plot.CIRCLE);	
			plot.show();
		}	
		
	}

	private void calculate_matches()
	{	
		float[] histpix = (float[])histImgSt.getPixels(1);	
		float[][] pdf_pixels = new float[14][histpix.length];
		
		
		double match_Poisson   = 0.0;
		double match_PoissonF  = 0.0;
		double match_Gauss     = 0.0;
		double match_Gamma     = 0.0;
		double match_Exp       = 0.0;
		double match_GaussF    = 0.0;
		double match_GammaF    = 0.0;
		double match_invGamma  = 0.0;
		double match_invBeta   = 0.0;
		double match_invGammaF = 0.0;
		double match_invBetaF  = 0.0;
		double match_Rubber    = 0.0;
		double match_lognorm   = 0.0;
		double match_lognormF  = 0.0;
		
		for(int y = 0; y < histheight; ++y) 
		{
			
			double delta_Poisson   = 0.0;
			double delta_PoissonF  = 0.0;
			double delta_Gauss     = 0.0;
			double delta_Gamma     = 0.0;
			double delta_Exp       = 0.0;
			double delta_GaussF    = 0.0;
			double delta_GammaF    = 0.0;
			double delta_invGamma  = 0.0;
			double delta_invBeta   = 0.0;
			double delta_invGammaF = 0.0;
			double delta_invBetaF  = 0.0;
			double delta_Rubber    = 0.0;
			double delta_lognorm   = 0.0;
			double delta_lognormF  = 0.0;
			
			double area_Poisson    = 0.0;
			double area_PoissonF   = 0.0;
			double area_Gauss      = 0.0;
			double area_Gamma      = 0.0;
			double area_Exp        = 0.0;
			double area_GaussF     = 0.0;
			double area_GammaF     = 0.0;
			double area_invGamma   = 0.0;
			double area_invBeta    = 0.0;
			double area_invGammaF  = 0.0;
			double area_invBetaF   = 0.0;
			double area_Rubber     = 0.0;
			double area_lognorm    = 0.0;
			double area_lognormF   = 0.0;
			
			double avg = avgs[y];
			double std = stdFitter.get_val(avg);//intercept1 + slope1 * avg;
			double variance = Math.pow(std,2);
			//Poisson Distribution
			double x0P = variance-avg;
			//Gamma Distribution
			double mG = avg+x0GM;
			double theta = variance/(mG);
			double ln_theta = Math.log(theta);
			double k = Math.pow(mG,2)/variance;
			double ln_Gk = lnGamma(k);
			//inverse Gamma Distribution
			double miG = avg+x0iGM;
			double alpha = +3.0+Math.sqrt(9.0-4.0*(2.0-miG/variance));
			double beta = miG*(alpha-1.0);
			double iGN = Math.pow(beta,alpha)*Math.exp(-lnGamma(alpha));
			//inverse Beta Distribution
			double miB = avg + x0iBM;
			double miB2 = Math.pow(miB,2);
			double miB3 = Math.pow(miB,3);
			//double iBalpha = (4*(miB3-miB2))/variance+miB;
			//double iBbeta = 1 + iBalpha/miB;
			double iBbeta = 2.0 + miB*(miB-1.0)/variance;
			double iBalpha = miB*(iBbeta-1.0);
			double iBN = lnEulerBeta(iBalpha,iBbeta);
			//lognormal Distribution
			double LNavg = avg + x0LN;
			double LNvar = Math.log(1+variance*Math.pow(LNavg,-2));
			double LNmu = Math.log(LNavg)-0.5*LNvar;
			double LNsig = Math.pow(LNvar,0.5);
			
			
			for(int x = 0; x < histwidth; ++x)
			{
				
				int xp = x + min_val;// + data_offset;
				double w = (double)histpix[x + y * histwidth];
				double[] ws = new double[14];
				
				double   w_Poisson   = probPoisson((int)(xp+x0Ps[y]+0.5),vars[y]);
				ws[0]  = w_Poisson;
				double w_PoissonF  = probPoisson((int)(xp+x0P+0.5),variance);
				ws[1]  = w_PoissonF;
				double w_Gauss     = probGauss(xp,avgs[y],stds[y]);
				ws[2]  = w_Gauss;
				double w_Gamma     = prob2Gamma(xp+x0GM, ks[y], thetas[y], lnThetas[y], lnGks[y]);
				ws[3]  = w_Gamma;
				double w_Exp       = probExp(xp,avgs[y]);
				ws[4]  = w_Exp;
				double w_GaussF    = probGauss(xp,avg,std);
				ws[5]  = w_GaussF;
				double w_GammaF    = prob2Gamma(xp+x0GM, k, theta, ln_theta, ln_Gk);
				ws[6]  = w_GammaF;
				double w_invGamma  = probiGamma(xp+x0iGs[y], alphas[y], betas[y], iGNs[y]);
				ws[7]  = w_invGamma;
				double w_invBeta   = probiBeta(xp+x0iBs[y], iBalphas[y], iBbetas[y], iBNs[y]);
				ws[8]  = w_invBeta;
				double w_invGammaF = probiGamma(xp+x0iGM, alpha, beta, iGN);
				ws[9]  = w_invGammaF;
				double w_invBetaF  = probiBeta(xp+x0iBM, iBalpha, iBbeta, iBN);
				ws[10]  = w_invBetaF;
				double w_Rubber = probRubber(xp,avgs[y],stds[y]);
				ws[11] = w_Rubber;
				double w_lognorm = probLN(xp+x0LNs[y],LNmus[y],LNvars[y]);
				ws[12] = w_lognorm;
				double w_lognormF = probLN(xp+x0LN,LNmu,LNvar);
				ws[13] = w_lognormF;
				
				for(int i=0; i < ws.length; ++i )
				{
					pdf_pixels[i][x + y * histwidth] = (float)ws[i];
					//transform to ws*ln(ws) for calculating pvaloffsets
					/*
					if(w==0)
					{	pvoffset += 1;}
					else
					{	pvoffset += w * Math.log(w);}
					
					if(ws[i]==0.0)
					{	pvoffsets[i] += 1.0;}
					else
					{	pvoffsets[i] += ws[i] * Math.log(ws[i]);}
					*/ 
				}
				
/*				
				pdf_pixels[ 0][x + y * histwidth] = (float)w_Poisson;
				pdf_pixels[ 1][x + y * histwidth] = (float)w_PoissonF;
				pdf_pixels[ 2][x + y * histwidth] = (float)w_Gauss;
				pdf_pixels[ 3][x + y * histwidth] = (float)w_Gamma;
				pdf_pixels[ 4][x + y * histwidth] = (float)w_Exp;
				pdf_pixels[ 5][x + y * histwidth] = (float)w_GaussF;
				pdf_pixels[ 6][x + y * histwidth] = (float)w_GammaF;
				pdf_pixels[ 7][x + y * histwidth] = (float)w_invGamma;
				pdf_pixels[ 8][x + y * histwidth] = (float)w_invBeta;
				pdf_pixels[ 9][x + y * histwidth] = (float)w_invGammaF;
				pdf_pixels[10][x + y * histwidth] = (float)w_invBetaF;
*/			
				area_Poisson   += w_Poisson;
				area_PoissonF  += w_PoissonF;
				area_Gauss     += w_Gauss;
				area_Gamma     += w_Gamma;
				area_Exp       += w_Exp;
				area_GaussF    += w_GaussF;
				area_GammaF    += w_GammaF;
				area_invGamma  += w_invGamma;
				area_invBeta   += w_invBeta;
				area_invGammaF += w_invGammaF;
				area_invBetaF  += w_invBetaF;
				area_Rubber    += w_Rubber;
				area_lognorm   += w_lognorm;
				area_lognormF  += w_lognormF; 
				
				
				delta_Poisson   += Math.abs(w_Poisson   - w);
				delta_PoissonF  += Math.abs(w_PoissonF  - w);
				delta_Gauss     += Math.abs(w_Gauss     - w);
				delta_Gamma     += Math.abs(w_Gamma     - w);
				delta_Exp       += Math.abs(w_Exp       - w);
				delta_GaussF    += Math.abs(w_GaussF    - w);
				delta_GammaF    += Math.abs(w_GammaF    - w);
				delta_invGamma  += Math.abs(w_invGamma  - w);
				delta_invBeta   += Math.abs(w_invBeta   - w);
				delta_invGammaF += Math.abs(w_invGammaF - w);
				delta_invBetaF  += Math.abs(w_invBetaF  - w);
				delta_Rubber    += Math.abs(w_Rubber    - w);
				delta_lognorm   += Math.abs(w_lognorm   - w);
				delta_lognormF  += Math.abs(w_lognormF  - w);
				
			}
					  
			
			// 0 .. Poisson, 1 .. Gauss, 2 .. Gamma, 3 .. neg. Exp 
			matches[y][ 0] = 1.0 - 2*delta_Poisson   / (area_Poisson   + 1 );
			matches[y][ 1] = 1.0 - 2*delta_PoissonF  / (area_PoissonF  + 1 );
			matches[y][ 2] = 1.0 - 2*delta_Gauss     / (area_Gauss     + 1 );
			matches[y][ 3] = 1.0 - 2*delta_Gamma     / (area_Gamma     + 1 );
			matches[y][ 4] = 1.0 - 2*delta_Exp       / (area_Exp       + 1 );
			matches[y][ 5] = 1.0 - 2*delta_GaussF    / (area_GaussF    + 1 );
			matches[y][ 6] = 1.0 - 2*delta_GammaF    / (area_GammaF    + 1 );
			matches[y][ 7] = 1.0 - 2*delta_invGamma  / (area_invGamma  + 1 );
			matches[y][ 8] = 1.0 - 2*delta_invBeta   / (area_invBeta   + 1 );
			matches[y][ 9] = 1.0 - 2*delta_invGammaF / (area_invGammaF + 1 );
			matches[y][10] = 1.0 - 2*delta_invBetaF  / (area_invBetaF  + 1 );
			matches[y][11] = 1.0 - 2*delta_Rubber    / (area_Rubber    + 1 );
			matches[y][12] = 1.0 - 2*delta_lognorm   / (area_lognorm   + 1 );
			matches[y][13] = 1.0 - 2*delta_lognormF  / (area_lognormF  + 1 );
			
			match_Poisson   += matches[y][ 0];
			match_PoissonF  += matches[y][ 1];
			match_Gauss     += matches[y][ 2];
			match_Gamma     += matches[y][ 3];	
			match_Exp       += matches[y][ 4];
			match_GaussF    += matches[y][ 5];
			match_GammaF    += matches[y][ 6];
			match_invGamma  += matches[y][ 7];
			match_invBeta   += matches[y][ 8];
			match_invGammaF += matches[y][ 9];
			match_invBetaF  += matches[y][10];
			match_Rubber    += matches[y][11];
			match_lognorm   += matches[y][12];
			match_lognormF	+= matches[y][13];		
		}
		
		int[] ordering = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
		double[] keys = {match_Poisson, match_PoissonF, match_Gauss, match_Gamma, match_Exp, match_GaussF,
			 match_GammaF, match_invGamma, match_invBeta, match_invGammaF, match_invBetaF,match_Rubber,match_lognorm, match_lognormF}; 	
		String[] titles = {"Poisson ", "PoissonF ", "Gauss ", "Gamma ", "Exp ", "GaussF ",
			 "GammaF ", "invGamma ","invBeta ", "invGammaF ","invBetaF ", "Rubber ", "LogNorm ", "LogNormF "};
		sort_by_first(keys,ordering,titles);
		if(create_PDF_maps)
		{
			ImagePlus pdfImg = NewImage.createFloatImage(WindowManager.makeUniqueName("pdf_maps.tif"),histwidth, histheight, ordering.length+1, NewImage.FILL_BLACK);
			ImageStack pdfImgSt = pdfImg.getStack();
			
			pdfImgSt.setPixels(histpix ,1);
			pdfImgSt.setSliceLabel("Data", 1);
			for(int i = 0; i < ordering.length; ++i)
			{
				pdfImgSt.setPixels(pdf_pixels[ordering[i]], i+2);
				pdfImgSt.setSliceLabel(titles[i] + keys[i], i+2);
			}
			pdfImg.setStack(pdfImgSt);
			if(!silent_mode) {	pdfImg.show();}
		}
	}
	
	private void make_prob_map()
	{
		
		if(create_table)
		{	table2 = new ResultsTable( );}
		else
		{	table2 = null;} //dont show zombies from last run
		
		
		model_min = avgs[ orderedYs[0] ]; //allways in the center of the hexagons
		
		if(modelBrightness < 0)
		{	
			double brightness_unit = (data_mean-model_min);
			model_min -= brightness_offset*brightness_unit; //barely sufficent for 8 ring
			brightness_unit = (data_mean-model_min);
			
			model_max = model_min - (modelBrightness * brightness_unit);	
			if(model_min < 0)
			{
				model_max -= model_min;
				model_min = 0;
			}		
		}
		else
		{	model_max = modelBrightness;}
		
		if( create_prob_map )
		{
			IJ.showStatus("generating Prob_Map.tif");
			pvaloffset = 0.0;
			float histpix[] = (float[])histImgSt.getPixels(1);
			
			for(int y = 0; y < histheight; ++y )
			{
				double sqr_sum = 0;
				for(int x = 0; x < histwidth; ++x)
				{	
					double fraction = histpix[x+y*histwidth];
					if(fraction > 0.0)
							//relative cases     * log probability
					sqr_sum += fraction * Math.log(fraction);
				}
				sqr_sum *= rel_Ws[y]; 
				pvaloffset += sqr_sum;
				
			}
			double pvfloor = pvalfloor ; //- pvaloffset ///TODO currently floor is fixed while offset is dynamic
			//IJ.log("pvaloffset: " + pvaloffset + "   pvfloor: " + pvfloor);
			probImg = NewImage.createFloatImage(WindowManager.makeUniqueName(probMapTitle),histwidth, modelGraylevels, 1, NewImage.FILL_BLACK);
			probImgSt = probImg.getStack();
			probImgSt.setSliceLabel(probTitle,1);
			float[] pixels = (float[])probImgSt.getPixels(1);
			double[] tmp_pixels = new double[pixels.length];;
			
			for(int mval = 0; mval < modelGraylevels; ++mval)
			{
				IJ.showProgress(mval,modelGraylevels);
				final double avg = model_min + (model_max-model_min)*mval/(modelGraylevels-1);
				double std = stdFitter.get_val(avg);//intercept1 + slope1 * avg;
				double variance = Math.pow(std,2);
				//Poisson Distribution
				double x0P = variance-avg;
				//Gamma Distribution
				double mG = avg+x0GM;
				double theta = variance/(mG);
				double ln_theta = Math.log(theta);
				double k = Math.pow(mG,2)/variance;
				double ln_Gk = lnGamma(k);
				//inverse Gamma Distribution
				double miG = avg+x0iGM;
				double alpha = 3.0+Math.sqrt(9.0-4.0*(2.0-miG/variance));
				double beta = miG*(alpha-1.0);
				double iGN = Math.pow(beta,alpha)*Math.exp(-lnGamma(alpha));
				//inverse Beta Distribution
				double miB = avg + x0iBM;
				double miB2 = Math.pow(miB,2);
				double miB3 = Math.pow(miB,3);
				//double iBalpha = (4*(miB3-miB2))/variance+miB;
				//double iBbeta = 1 + iBalpha/miB;
				double iBbeta = 2.0 + miB*(miB-1.0)/variance;
				double iBalpha = miB*(iBbeta-1.0);
				double iBN = lnEulerBeta(iBalpha,iBbeta); 
				
				//lognormal Distribution
				double LNavg = avg + x0LN;
				double LNvar = Math.log(1+variance*Math.pow(LNavg,-2));
				double LNmu = Math.log(LNavg)-0.5*LNvar;
				double LNsig = Math.pow(LNvar,0.5);
				
				
				if(table2 != null)
				{	
					table2.setValue("avg.", mval, avg);
					table2.setValue("std.", mval, std);
					table2.setValue("var.", mval, variance);
					table2.setValue("x0P", mval, x0P);
					table2.setValue("x0G", mval, x0GM);
					table2.setValue("k", mval, k);
					table2.setValue("theta", mval, theta);
					table2.setValue("x0iG", mval, x0iGM);
					table2.setValue("iGalpha", mval, alpha);
					table2.setValue("iGbeta", mval, beta);
					table2.setValue("x0iB", mval, x0iBM);
					table2.setValue("iBalpha", mval, iBalpha);
					table2.setValue("IBbeta", mval, iBbeta);	
					table2.setValue("x0LN", mval, x0LN);
					table2.setValue("LNmu", mval, LNmu);
					table2.setValue("LNsig", mval, LNsig);			
				}
				
				for(int i=0; i < histwidth; ++i)
				{
					int x = i + min_val;
					double xf0 = x;
					double xf = x; 
					double prob = 0.0;
					if( interpolate && (avg <= avgs[ orderedYs[ histheight-1 ] ]) && (avg >= avgs[ orderedYs[ 0 ] ]) )
					{
						prob = probInterpolate(x, avg, std);
					}
					else
					{
						switch(probChoice)
						{
							case 1:
								x = (int)(x+x0P+0.5);
								prob = probPoisson(x, variance);
							break;
							case 2:
								xf = (xf0+x0GM);
								prob = prob2Gamma(xf, k, theta, ln_theta, ln_Gk);
							break;
							case 3:
								prob = probGauss(x, avg, std);
							break;
							case 4:
								prob = probExp(x, avg);
							break;
							case 5:
								prob = probHisto(x, avg, std);
							break;
							case 6: 
								xf = (xf0+x0iGM);
								prob = probiGamma(xf,alpha,beta,iGN);
							break;
							case 7:
								xf = (xf0+x0iBM);
								prob = probiBeta(xf,iBalpha,iBbeta,iBN);
							break;
							case 8:
								prob = probRubber(x,avg,std);
							break;	
							case 9:
								xf = (xf0+x0LN);
								prob = probLN(xf,LNmu,LNvar);
							break;
							default:
								IJ.error("Unknown Choice of PDF for creating " + probImg.getTitle());
								return;
						}
					}
					tmp_pixels[i + mval * histwidth] = prob;
				}
			}
			
			for(int mval = 0; mval < modelGraylevels; ++mval)
			{
				final double avg = model_min + (model_max-model_min)*mval/(modelGraylevels-1);
				if( ( interpolate && (avg <= avgs[ orderedYs[ histheight-1 ] ]) ) || (probChoice == 5) || (probChoice == 8) ) //
				{
					double sum = 0.0;
					for(int i=0; i < histwidth; ++i)
					{
						sum += tmp_pixels[i + mval * histwidth];
					}
					final double inv_sum = 1.0/sum;
					for(int i=0; i < histwidth; ++i)
					{
						tmp_pixels[i + mval * histwidth] *= inv_sum; 
					}
				}
			}
			for(int ind=0; ind < pixels.length; ++ind)
			{
				pixels[ind] = (float)tmp_pixels[ind];
			}
			
			probImg.setStack(probImgSt);
			if(gb_sigma > 0.0)  
			{
				GaussianBlur gb = new GaussianBlur();
				FloatProcessor fp = (FloatProcessor)probImg.getProcessor();
				gb.blur1Direction(fp,gb_sigma,0.01,true,0);
			}
			
			for(int ind=0; ind < pixels.length; ++ind)
			{
				double prob = (double)pixels[ind];
				double logprob = Math.log(prob) - pvaloffset;
				if (logprob<pvaloffset) //a regularization.
				{	logprob=pvaloffset*(1.0d+Math.tanh((logprob-pvaloffset)/pvaloffset));}
				pixels[ind] = (float)logprob;
			}
			
			if(!silent_mode) 
			{	
				probImg.show();
				if(table2 != null)
				{	
					table2.showRowNumbers(true);
					table2.show("Fits " + (Idata!=null?Idata.getTitle():histImg.getTitle()));
				}
			}	
		}
	}
	
	private void setHexcenter(int q, int r)
	{
		//in principle negative centers are allowed
		cq = q;
		cr = r;
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
	
	private void make_hex_mask()
	{
		hex_mask = new boolean[4 * rH * rH];
		for(int ind = 0; ind < hex_mask.length ; ++ind)
		{	
			int q = ind%(2*rH);
			int r = ind/(2*rH);		
			hex_mask[ind] = iswithinHex(q,r,rH);		
		
		}
	}
	
	
	public boolean hexagonaltrap(int[] cube, int N)
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
	
	public boolean sort_by_first(double[] keys, int[] ordering, String[] titles)
	{
		/*
		IJ.log("sorting now");
		IJ.log("keys: " + keys[0] + " " + keys[1] + " " + keys[2] + " " 
							+ keys[3] + " " + keys[4] + " " + keys[5] );
		IJ.log("ordering: " +ordering[0] + " " + ordering[1] + " " + ordering[2] + " " 
							+ ordering[3] + " " + ordering[4] + " " + ordering[5] );					
		IJ.log("titles: " + titles[0] + " " + titles[1] + " " + titles[2] + " " 
							+ titles[3] + " " + titles[4] + " " + titles[5] );	
		*/
		boolean ordered = true;
		for(int i = 0; i < keys.length-1; ++i)
		{
			for( int j = i+1; j < keys.length; ++j)
			{
				if(keys[i] < keys[j])
				{//swap the elements in keys and ordering
					//IJ.log("swaping " + i + " " + j + " ordering: " + ordering[i] + " <-> " + ordering[j] + " titles: " + titles[i] + " <-> " + titles[j] );
					double tmp_key = keys[i];
					int tmp_ord = ordering[i];
					String tmp_title = titles[i];
					
					keys[i] = keys[j];
					ordering[i] = ordering[j];
					titles[i] = titles[j];
					
					keys[j] = tmp_key;
					ordering[j] = tmp_ord;
					titles[j] = tmp_title;
					ordered = false;
					/*
					IJ.log("keys: " + keys[0] + " " + keys[1] + " " + keys[2] + " " 
							+ keys[3] + " " + keys[4] + " " + keys[5] );
					IJ.log("ordering: " +ordering[0] + " " + ordering[1] + " " + ordering[2] + " " 
							+ ordering[3] + " " + ordering[4] + " " + ordering[5] );
					IJ.log("titles: " + titles[0] + " " + titles[1] + " " + titles[2] + " " 
							+ titles[3] + " " + titles[4] + " " + titles[5] );			
					*/	
				}
			}
		}
		//IJ.log("all sorted?");
		return ordered;
	}
	
	public boolean order_by_key(double[] keys, int[] ordering)
	{
		double[] mykeys = new double[keys.length];
		for(int i = 0; i < keys.length; ++i)
		{
			mykeys[i] = keys[i];
		}
		boolean ordered = true;
		for(int i = 0; i < keys.length-1; ++i)
		{
			for( int j = i+1; j < keys.length; ++j)
			{
				if(mykeys[i] > mykeys[j])
				{//swap the elements in keys and ordering
					double tmp_key = mykeys[i];
					int tmp_ord = ordering[i];
					
					mykeys[i] = mykeys[j];
					ordering[i] = ordering[j];
					
					mykeys[j] = tmp_key;
					ordering[j] = tmp_ord;
					ordered = false;
					
				}
			}
		}
		return ordered;
	}
	
		
	private boolean get_input()
	{
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog(called_by);
		gd.setSmartRecording(true);
		int[] idArray = WindowManager.getIDList(); // list of all opened images (IDs)
		int idlen = 0;
		if(idArray != null)
		{	idlen = idArray.length;}
		String[] titleArray = new String[idlen + 2]; // titles of opened images
		String[] histoArray = new String[idlen + 2];
		
		for (int i = 0; i < idlen; ++i)
		{	
			titleArray[i] = WindowManager.getImage(idArray[i]).getTitle();
			histoArray[i+1] = WindowManager.getImage(idArray[i]).getTitle();
		}
		titleArray[idlen] = "<open>";
		titleArray[idlen+1] = "<none>";
		histoArray[0] = "<create>";
		histoArray[idlen+1] = "<open>";
		if(dataTitle==null && WindowManager.getCurrentImage()!= null)	
		{	dataTitle = WindowManager.getCurrentImage().getTitle(); }
		gd.addChoice("Data input", titleArray, dataTitle);
		gd.addCheckbox("virtual smoothing", virtual_smoothing);
		
		gd.addNumericField("graphene C-C bondlength" , bondlength, 0);
		gd.addNumericField("modelnum (<1 none)", modelnum, 0);
		gd.addNumericField("Minimal width of Histogramm", min_histwidth,0);
		gd.addCheckbox("show equivalent sites", create_wedge);
		gd.addChoice("Hex_Histogramm", histoArray, histTitle);
		gd.addCheckbox("symmetrized Histogramm", sym_Histo);
		gd.addCheckbox("show PDF maps", create_PDF_maps);
		gd.addCheckbox("show Tables ", create_table);
		gd.addCheckbox("show Plot", plot_std_vs_avg);
		gd.addMessage("0 .. none, 1 .. basic, 2 .. statistics, 3 .. fits, 4 .. debug");
		gd.addNumericField("output_level" , output_level, 0);
		//gd.addNumericField("data_offset" , data_offset, 0);
		gd.addChoice("fit model", stdFitter.get_all_models(), fitname);
		gd.addChoice("Log_Prop_map", probArray, probTitle);
		gd.addCheckbox("interpolate propabilities", interpolate);
		gd.addCheckbox("keep step", keep_step);
		gd.addNumericField("threshold for histogramm: ", minfactor, 0);
		gd.addNumericField("Model gray levels: ", modelGraylevels, 0);
		gd.addNumericField("Gaussian Blur sigma (<0 none): ", gb_sigma, 3);
		gd.addNumericField("Equivalent max data in model1(x<0 x*mean)", modelBrightness, 3);
		
		gd.addNumericField("regularization for log probs",pvalfloor, 3);
		
		gd.addCheckbox("quit after operation",!interactive_mode);
		if(!silent_mode)
		{
			gd.showDialog();
			if( gd.wasCanceled() )
			{	return false;	}
		}
		dataTitle = gd.getNextChoice();
		virtual_smoothing = gd.getNextBoolean();
		bondlength = (int) gd.getNextNumber();
		modelnum = (int) gd.getNextNumber();
		min_histwidth = (int) gd.getNextNumber();
		create_wedge = gd.getNextBoolean();
		histTitle = gd.getNextChoice();
		sym_Histo = gd.getNextBoolean();
		create_PDF_maps = gd.getNextBoolean();
		create_table = gd.getNextBoolean();
		plot_std_vs_avg = gd.getNextBoolean();
		output_level = (int) gd.getNextNumber();
		min_val = 0; //data_offset;
		fitname = gd.getNextChoice();
		probChoice = (int) gd.getNextChoiceIndex();
		probTitle = probArray[probChoice];
		create_prob_map = (probChoice != 0);
		interpolate = gd.getNextBoolean();
		keep_step = gd.getNextBoolean();
		minfactor = (int) gd.getNextNumber();
		modelGraylevels = (int) gd.getNextNumber();
		gb_sigma = gd.getNextNumber();
		modelBrightness = gd.getNextNumber();
		pvalfloor = gd.getNextNumber();
		interactive_mode = !gd.getNextBoolean();
		
		return true;
	}
	
	private boolean validate_input()
	{
		//we dont want to default to prior runs
		Idata = null;
		IdataSt = null;
		histImg = null;
		histImgSt = null;
		
		boolean status = true;
		if( !histTitle.equals("<create>") )
		{
			if( !dataTitle.equals("<none>") )
			{
				if(output_level > 3)
				{
					IJ.log("Data input must be \"<none>\" if using an existing Hex_Histgramm");
					return false;
				}
				else
				{
					dataTitle = "<none>";
				}
			}
			if(histTitle.equals("<open>"))
			{
				histTitle = IJ.getFilePath("open Histogramm");
				if(histTitle == null)
				{	return false;}
				histImg = new ImagePlus(histTitle);
				histImg.show();
				histTitle = histImg.getTitle();				
			}
			else if (! (histTitle.startsWith("<") && histTitle.endsWith(">")))
			{	histImg = WindowManager.getImage(histTitle);}
			if( histImg != null)
			{	histImgSt = histImg.getStack();}
			
			if( (histImg != null) && (histImg.getType() != ImagePlus.GRAY32) )
			{
				status = false;
				IJ.log("histogramm data is not GRAY32, please convert " + histTitle);
			}
			if( create_wedge )
			{
				if(output_level > 3)
				{	IJ.log("WARNING: map of equivalent pixels is only available when creating a new Hex_Histogramm");}
				create_wedge = false;
			}
			
		}
		
		if(dataTitle.equals("<open>"))
		{
			dataTitle = IJ.getFilePath("open Idata");
			if(dataTitle == null)
			{	return false;}
			Idata = new ImagePlus(dataTitle);
			Idata.show();
			dataTitle = Idata.getTitle();				
		}
		else if (! (dataTitle.startsWith("<") && dataTitle.endsWith(">")))
		{	Idata = WindowManager.getImage(dataTitle);}
		if( (Idata != null) && (Idata.getType() != ImagePlus.GRAY16) )
		{
			IJ.log("raw data is not GRAY16, please convert " + dataTitle);
			status = false;
		}
		if( Idata != null) //check the stuff that is relevant for processing the raw data
		{
			IdataSt = Idata.getStack();
			modelsize = Idata.getWidth();
			modelarea = modelsize * modelsize;
			impDepth = Idata.getStackSize();
		
			if(modelsize != Idata.getHeight())
			{
				IJ.log(dataTitle + " must be square shaped");
				status = false;
			}
			
			if( (bondlength < 2) || (bondlength % 2 == 1) )
			{	
				status = false;
				IJ.log("C-C bondlength must be positive and even: " + bondlength);
			}
			
			if( (modelsize < 1) || ( modelsize % (2*bondlength) != 0 ) )
			{	
				status = false;
				IJ.log("modelsize must be a positive even mutliple of C-C bondlength: " + modelsize);
			}
			
			rH = modelsize/2;
			setHexcenter( modelsize/2, modelsize/2 );
			make_hex_mask();	
		}
		
		//if(data_offset < 0)
		//{	IJ.log("WARNING: data offset should typically be >= 0 to account for an earlier background subtraction");}
		
		if(modelGraylevels < 2)
		{
			IJ.log("There cannot be less than 2 grayvalus in a model");
			status = false;
			modelGraylevels = 2;
		}
		else if(modelGraylevels < 8)
		{
			IJ.log("WARNING: Less than 8 graylevels are really not recommended");
		}
		
		if(minfactor < 0)
		{
			IJ.log("Negative histogram threshold values dont make sense as probabilities");
			status = false;
			minfactor = 1;
		}
		else if (minfactor == 0)
		{
			IJ.log("WARNING: threshhold is 0. There may be plenty of NaNs and Infs ahead");
		}
		
		
		return status;
	}
	
	private void log_results()
	{
		if(output_level > 0)
		{
			if( Idata != null)
			{
				IJ.log("Input: " + Idata.getTitle() + "  min: " + real_min_val + "@" + min_pos 
													+ "  max: " + max_val + "@" + max_pos
													+ "  mean: "+ data_mean );
				IJ.log("modelsize: " + modelsize + "  bondlength: " + bondlength);
			}
			if(	modelnum > 0)
			{
				IJ.log("created " + modelnum + " models");
				IJ.log("model min,max " + model_min + "," + model_max + "  Graylevels: " + modelGraylevels);
			}
			
			if(histTitle.equals("<create>"))
			{	
				IJ.log("created Histogramm" + histImg.getTitle()); 
				if(output_level > 3)
				{
					for(int y = 0; y < histheight; ++y)
					{	IJ.log("linesum" + y + ": " + linesums[y] + "  fraction: " + rel_Ws[y] );}
				}	
			}
			else 
			{	IJ.log("loaded Histogramm " + histImg.getTitle());}
			if(create_prob_map)
			{
				IJ.log("created " + probImg.getTitle() + " from " + probTitle + " pvaloffset: " + pvaloffset + " blur: " + gb_sigma);
				IJ.log("interpolation: " +  interpolate + "   keep_step: " + keep_step  + "   minfactor: " + minfactor);
				IJ.log("ModelGraylevels(vertical) " + modelGraylevels + " between " + model_min + " and " + model_max);
				IJ.log("data values(horizontal) from " + (min_val) + " to " + (max_val) );
			}
		
		}
		if(output_level > 1)
		{
			IJ.log("Weighted offsets: Poisson: " + x0PM + 
						            " Gamma: " + x0GM + 
						            " invGamma: " + x0iGM +
						            " invBeta: " + x0iBM);
			
			//IJ.log("offset:" + data_offset);
		}
		
		if(output_level > 2)
		{
			IJ.log("Regression std vs. avg: " + stdFitter.get_model() + "  rms: " + stdFitter.get_merit() );
			// 0 .. Poisson, 1 .. Gauss, 2 .. Gamma, 3 .. Exp, 4 .. GaussF, 5 .. GammaF
			if(output_level > 3)
			{
				IJ.log("match is defined as 1-2*(common area)/(total area)");
			}
		}
		
		return;	
	}
		
	public ImagePlus get_prob_map(String dT, String PDF, String pmT, boolean symH, int outlvl, int bl, int mgl,
	 int mf, double mBs, double pvflr, double gbs, boolean intpol, boolean keepst)
	{
			silent_mode = true;
			interactive_mode = false;
			create_prob_map = true;
			create_table = false;
			virtual_smoothing = false;
			interpolate = intpol;
			keep_step = keepst;
			output_level = outlvl;
			modelnum = 0;
			dataTitle = dT;
			bondlength = bl;
			probTitle = PDF;
			probMapTitle = pmT;
			histTitle = "<create>";
			sym_Histo = symH;
			output_level = 0;
			minfactor = mf;
			modelGraylevels = mgl;
			modelBrightness = mBs;
			pvalfloor = pvflr;
			gb_sigma = gbs;
			run("remote");
			return probImg;
	}
	
	
	
	public void run(String arg)
	{
		called_by = arg;
		boolean redo = false;
		do
		{
			redo = false;
			if( !get_input() )
			{	break;	}
			if( !validate_input() )
			{	
				redo = true;
				silent_mode = false;
				continue;	
			}
			make_unique_pixels();
			if(create_wedge)
			{	show_inequivalent_pixels();}
			setup_histImg();
			if(!dataTitle.equals("<none>"))
			{	
				populate_histImg();
			}
			if( (output_level > 1) || create_PDF_maps || create_prob_map || create_table || plot_std_vs_avg || (modelnum > 0) )
			{	do_statistics();}
			
			if( (output_level > 2) || create_PDF_maps )
			{	calculate_matches();}
			if( create_prob_map || (modelnum > 0))
			{	make_prob_map();}
			if(modelnum > 0 && (!silent_mode) )
			{	
				make_models("Models.tif", modelnum);
				make_asym_models("Asym_Models.tif", modelnum);	
			}
			log_results();
		}while( interactive_mode || redo);			
		return;
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
	
	private double lnEulerBeta(double a, double b)
	{
		return lnGamma(a) + lnGamma(b) - lnGamma(a+b);
	}
	     
	private double probPoisson(int x, double avg)
	{
		if(x < 0) {	return 0;}
		if( avg <= 0.0)
		{	return (x==0)?1.0:0.0;}
		if( x == 0)
		{	return Math.exp(-avg);}
		int kk = x;
		double f0 = Math.log(avg);
		double f = (double)x*f0 - avg;
		while (kk > 1)	{	f -= Math.log( (double)(kk--) );}
		return Math.exp(f);
	}

	private double probExp(int x, double avg)
	{
		if(avg <= 0.0)
		{	return (x==0)?1.0:0.0;}
		if(x < 0)
		{	return 0.0;}
		return Math.exp( -x/avg - Math.log(avg) );
	}
	
	private double probInterpolate(int x, double avg, double std) //no std needed here
	{
		x -= min_val;
		float histpix[] = (float[])histImgSt.getPixels(1);
		double result = 0;
		int lower_y = 0;
		int upper_y = histheight-1;
		for(int y = 0; y < histheight; ++y)
		{
			
			if( (avgs[y] <= avg) && (avgs[lower_y] < avgs[y]) )
			{	lower_y = y;}
			int yr = histheight-1-y; //reverse
			if( (avgs[yr] >= avg) && (avgs[upper_y] > avgs[yr])  )
			{	upper_y = yr;}
		}
		double lk = stds[lower_y]/std;
		double lkk = (x > avg)?1/lk:lk;
		double minval_ly = minfactor * lkk * Math.pow(linesums[lower_y],-1);
		
		if(lower_y == upper_y) //exact match
		{
			result = Math.max(histpix[ x + lower_y*histwidth],minval_ly);
		}
		else //interpolate lower and upper
		{
			double uk = stds[upper_y]/std;
			double ukk = (x > avg)?1/uk:uk;
			double minval_uy = minfactor * ukk * Math.pow(linesums[upper_y],-1);
			
			double pt_lower = Math.max(histpix[ x + lower_y*histwidth],minval_ly);
			double pt_upper = Math.max(histpix[ x + upper_y*histwidth],minval_uy);
			double mixing = (avg-avgs[lower_y])/(avgs[upper_y]-avgs[lower_y]);
			result = ( (1.0-mixing)*pt_lower + (mixing)*pt_upper );	
		}
		return	result;
	}
	
	private double probRubber(int x, double avg, double std)
	{
		float histpix[] = (float[])histImgSt.getPixels(1);
		double area = 0;
		int ptsy = 0;
		for(int y = 0; y < histheight; ++y)
		{
			
			double k = stds[y]/std;
			double k0 = k;
			if(keep_step && (x < modes[y]) )
			{ 
				double filling = Math.sqrt(histpix[x-min_val + y * histwidth] / peak_vals[y]);
				k = 1 + filling * (k - 1 );
				
			}	
			//kk causes falling vertical cut on right edge of Prob_Map
			//and raising vertical cut on left edge
			double kk = (x > avg) ? 1/k0 : k0;
			double d = avgs[y] - k*avg;
			double minval_y = minfactor*kk*Math.pow(linesums[y],-1);
			int m1 = (int)Math.floor( k*(x-0.5) + d - min_val );
			int m2 = (int)Math.ceil(  k*(x+0.5) + d - min_val );
			if(m1 < 0) {m1 = 0;}
			else if(m1 >= histwidth) {m1 = histwidth-1;}
			if(m2 < 0) {m2 = 0;}
			else if(m2 >= histwidth) {m2 = histwidth-1;}
			int mpts = 0;
			double ay = 0.0;
			for(int m = m1; m <=m2; ++m)
			{
				ay +=  Math.max(k*histpix[m + y*histwidth] , minval_y); //k*
				++mpts;//+=1.0/k;
			}
			
			area += ( (mpts > 0) ? ay/mpts : minval_y );
			++ptsy;	
		}
		return area/ptsy;
	}
	
	
	private double probHisto(int x, double avg, double std/*, double avg_l, double avg_h*/)
	{
		x -= min_val;
		float histpix[] = (float[])histImgSt.getPixels(1);
		double area = 0;
		
		double[] prox_Ws = new double[histheight];
		double prox_norm = 0.0;
		for(int y = 0; y < histheight; ++y)
		{
			prox_Ws[y] = Math.pow(0.001 + Math.abs(avg - avgs[y]) , -1);
			prox_norm += prox_Ws[y];	
		}
		for(int y = 0; y < histheight; ++y)
		{
			prox_Ws[y]/=prox_norm;
		}
			
		for(int y = histheight-1; y < histheight; ++y)
		{
			double minval_y = Math.pow(linesums[y],-1);
			double shift_avg = avg - avgs[y];
			double shG = shift_avg * Math.sqrt(2/Math.PI);
			double var_y = (int)(std*std - vars[y]);
			double std_y = Math.sqrt(Math.max(0.25,var_y));
			double area_y = 0;
			double pts_y = 0;
			if( (shift_avg >= 0) /*&& (var_y > -1)*/ )
			{
				double sa = shift_avg - 0.5;
				// Some Gamma correction code
				// + Math.pow(shift_avg*(avgs[histheight-1]-avgs[0])/(model_max-model_min),1.5 );
				for(int m = 0; m < histwidth; ++m)
				{
					final int n = x - m;
					if( (m >= 0) && (n>0) )
					{
						//final double pG = probGauss(n , 0 , shG ); 	 //not so smooth continuation	
						final double pP = probExp(n , sa/*, std_y*/ ); //quite smooth continuation
						//final double pP = probGauss(n , sa, std_y ); //kink as compared to interpolation
						//final double pP = probPoisson(n , shift_avg );//kink in PropMap
						//double pP = probGamma(n , shift_avg, std_y ); //major discontiuity
						area_y += pP*Math.max(histpix[m + y*histwidth],minval_y);
						pts_y += pP;
					}
				}
			}
			if(pts_y > 0)
			{
				area_y /= pts_y;
			}
			else
			{
				area_y = minval_y;
			}	
			area += area_y * prox_Ws[y];
		}
		return area;
	}

	private double probGauss(int x, double avg, double std)
	{
		return 	Math.exp( -0.5*Math.pow( ((double)x-avg)/std , 2 )  ) * 
				Math.pow( 2*Math.PI, -0.5 ) / std;
	}

	private double probiGamma(double x, double alpha, double beta, double igN)
	{
		if(x < 0) {return 0;}
		return igN*Math.pow(x,-(alpha+1))*Math.exp(-beta/x);
	}

	private double probiBeta(double x, double alpha, double beta, double ibN)
	{
		if(x <= 0) {return 0.0;}
		return Math.exp( (alpha-1.0)*Math.log(x) - (alpha+beta)*Math.log(x+1.0) - ibN );
	}

	private double prob2Gamma(double x, double k, double theta, double ln_theta,double ln_Gk)
	{	
		if(x < 0)
		{	return 0;}
		if(theta <= 0.0)
		{	return (x==0)?1.0:0.0; } 
		double g1 = 1.0;
		if( (x==0) && (k <= 1.0) ) //just calculate the missing area from discrete distribution
		{
			double dg1 = 0.0;
			double x0 = 1;
			do
			{
				dg1 = prob2Gamma(x0, k, theta, ln_theta, ln_Gk);
				g1-= dg1;
				++x0;
			}while( (dg1 >= 0.000000001) && ( x0 < (k-1.0)*theta) );
			return g1;
		}
		return Math.exp(	(k-1.0) * Math.log( x )
						- k * ln_theta
						- ln_Gk - x / theta
					   );
	}
	
	private double probGamma(int x, double avg, double std)
	{	
		double variance = Math.pow(std,2);
		double theta = variance/avg;
		double k = Math.pow(avg,2)/variance;
		double ln_Gk = lnGamma(k);	
		if(theta <= 0.0)
		{	return (x==0)?1.0:0.0; } 
		double g1 = 1.0;
		if( (x==0) && (k <= 1.0) ) //just calculate the missing area from discrete distribution
		{
			double dg1 = 0.0;
			int x0 = 1;
			do
			{
				dg1 = probGamma(x0, avg, std);
				g1-= dg1;
				++x0;
			}while( (dg1 >= 0.000000001) && ( (double)x0 < (k-1.0)*theta) );
			return g1;
		}
		return Math.exp(	(k-1.0) * Math.log( (double)x )
						-k    * Math.log(theta)
						-ln_Gk - ( (double)x ) / theta
					   );
	}
	
	private double probLN(double x, double mu, double sig2)
	{
		if(x <= 0.0)	return 0.0;
		return Math.pow(x*Math.pow(2*Math.PI*sig2,0.5),-1)*Math.exp(-Math.pow(Math.log(x)-mu,2)/(2*sig2));
	}

}
