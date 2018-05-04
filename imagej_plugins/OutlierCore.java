import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import java.awt.event.*;
import ij.plugin.filter.*;
import ij.gui.DialogListener;
import ij.plugin.filter.PlugInFilterRunner;
import java.util.Arrays;

class OutlierCore implements DialogListener
{
	int radius1 = 2;
	int radius2 = 1;
	int weight = 1;
	int threshold = 50;
	int passes = 1;
	int mode = 0;
	final String[] modes = {"median","mean","min","max"};
	boolean hex_pixels = true;
	boolean filter_bright = true;
	boolean filter_dark   = true;
	boolean create_report = false;	 
	int cq = 0;
	int cr = 0;
	int cx = 0;
	int cy = 0;
	int cz = 0;
	int crH = 0;
	int imp_width = 0;
	int imp_height = 0;
	int flags;
	int DONE = 4096;
	
	OutlierCore(int flg) //dont call that with flg == DONE == 4096
	{
		flags = flg;
	};
	
	public int showDialog(final ImagePlus imp)
	{
		return showDialog(imp, "Removing Outliers", null );
	}
	
	public int showDialog(final ImagePlus imp, final String command)
	{
		return showDialog(imp, command, null);
	}
	
	public int showDialog(final ImagePlus imp, final String command, final PlugInFilterRunner pfr)
	{
		setimg(imp);
		GenericDialog gd = ((pfr != null)? new GenericDialog( command ) : new NonBlockingGenericDialog( command )); 
		gd.addMessage( imp.getTitle() );
		gd.addCheckbox("hex kernel & pixels", hex_pixels);
		//gd.addCheckbox("hex domain", hex_domain);
		gd.addNumericField("radius (detection)", radius1, 0);
		gd.addNumericField("radius (replacement)", radius2, 0);
		gd.addNumericField("weight of central pixel",weight,0);
		gd.addMessage("use threshold = 0 for smoothing");
		gd.addNumericField("threshold", threshold, 0);
		gd.addCheckbox("filter too bright", filter_bright);
		gd.addCheckbox("filter too dark", filter_dark);
		gd.addNumericField("passes", passes, 0);
		gd.addChoice("filter mode", modes, modes[mode] );              
		if(pfr != null)
		{	
			gd.addPreviewCheckbox(pfr);
			gd.addDialogListener(this);	
		}
		else
		{	gd.addMessage("Preview is available in Outlier_Remover");}
		gd.addCheckbox("create report", create_report);
		
		gd.showDialog();
        if (gd.wasCanceled())
		{	return DONE;}
		hex_pixels = gd.getNextBoolean();
		//hex_domain = gd.getNextBoolean();
		radius1 = (int)gd.getNextNumber();
		radius2 = (int)gd.getNextNumber();
		weight = (int)gd.getNextNumber();
		threshold = (int)gd.getNextNumber();
		filter_bright = gd.getNextBoolean(); 
		filter_dark = gd.getNextBoolean(); 
		passes = (int)gd.getNextNumber();
		mode = gd.getNextChoiceIndex();
		hex_pixels = gd.getNextBoolean();
		IJ.register(getClass());
		if( modes[mode].equals("min") )
		{
			filter_bright = true;
			filter_dark = false;
		}
		else if( modes[mode].equals("max") )
		{
			filter_bright = false;
			filter_dark = true;
		}
		if( (!filter_bright) && (!filter_dark) && (threshold > 0) )
		{	return DONE;}
		if ( ( !gd.invalidNumber() ) && (radius1 > 0) && (radius2 > 0) && (passes > 0) && (threshold >= 0) && (weight >= 0) )
		{	
			if(create_report)
			{	do_report(command, imp);}
			return flags;
		}
		else
		{	
			return DONE;
		}

	}
	
	public boolean dialogItemChanged(final GenericDialog gd, final AWTEvent e) 
	{
		hex_pixels = gd.getNextBoolean();
		//hex_domain = gd.getNextBoolean();
		radius1 = (int)gd.getNextNumber();
		radius2 = (int)gd.getNextNumber();
		weight = (int)gd.getNextNumber();
		threshold = (int)gd.getNextNumber();
		filter_bright = gd.getNextBoolean(); 
		filter_dark = gd.getNextBoolean(); 
		passes = (int)gd.getNextNumber();
		mode = gd.getNextChoiceIndex();
		if( modes[mode].equals("min") )
		{
			filter_bright = true;
			filter_dark = false;
		}
		else if( modes[mode].equals("max") )
		{
			filter_bright = false;
			filter_dark = true;
		}
		
		Checkbox brightCheckbox = (Checkbox)gd.getCheckboxes().get(1);
		brightCheckbox.setState(filter_bright);
		Checkbox darkCheckbox = (Checkbox)gd.getCheckboxes().get(2);
		darkCheckbox.setState(filter_dark);
		
		if((!filter_bright) && (!filter_dark) && (threshold > 0))
		{	return false;}
        return ( !gd.invalidNumber() ) && (radius1 > 0) && (radius2 > 0) && (passes > 0) && (threshold >= 0) && (weight >= 0);
    }
	
	private void do_report(final String command, final ImagePlus imp)
	{
		IJ.log(command);
		IJ.log(imp.getTitle());
		IJ.log(hex_pixels ? "hex kernel & pixels" : "square kernel & pixels" );
		IJ.log("radius(detection) : " + radius1 );
		IJ.log("radius(replacment) : " + radius2 );
		IJ.log("central weight : " + weight );
		IJ.log("threshold: " + threshold );
		IJ.log("filter too bright: " + filter_bright);
		IJ.log("filter too dark: " + filter_dark);
		IJ.log("mode: " + modes[mode]);
		return;
	}
	
	
	
	public void setup(final int r1, final int r2, final int wght,
				final int thold, final int pses,
				final int md, final boolean hp,
				final boolean fb, final boolean fd/*,
				final int q, final int r,
				final int x, final int y, final int z,
				final int rH, final int iw, final int ih*/ )
	{
		radius1 = r1;	radius2 = r2; weight = wght;
		threshold = thold;	passes = pses;
		mode = md;	hex_pixels = hp;
		filter_bright = fb;		filter_dark = fd;
		/*
		cq = q;	cr = r;
		cx = x;	cy = y;	cz = z;
		crH = rH;	imp_width = iw;	imp_height = ih;
		*/
		//IJ.log("OutlierCore received defaults");
		return;
	}
	
	public int setimg( final ImagePlus imp )
	{
		if(imp == null)
		{	
			IJ.noImage();
			return DONE;
		}
		else
		{
			imp_width = imp.getWidth();
			imp_height = imp.getHeight();
			int rH = imp_width/2;
			if(imp_height/2 < rH)
			{	rH = imp_height/2;}
			setHexcenter(imp_width/2,imp_height/2,rH);
			//IJ.log("OutlierCore initialized by image");
			return flags;
		}
	}
	
	private void setHexcenter(final int q, final int r, final int rr)
	{
		//in principle negative centers are allowed
		cq = q;
		cr = r;
		cx = q - ( r - (r&1) ) / 2;
		cz = r;
		cy = -cx - cz;
		crH = rr;
	}
	
	
	public void run( short[] pixels) 
	{
		final int th = threshold;
		final int r1 = radius1;
		final int r2 = radius2;
		final boolean reuse = ( r1==r2 );
		final boolean squares = (!hex_pixels);
		//short[] pixels = (short[])ip.getPixels();
        
        
        Filter filter = null;
        if(hex_pixels) 
        {	filter = new HexFilter();}
        else
        {	filter = new SquareFilter();}
        for(int i = 0; i < passes; ++i)
        {
			filter.setOldPix( pixels );
			for(int ind = 0; ind < pixels.length; ++ind)
			{
				final int q = ind % imp_width;
				final int r = ind / imp_width;
				if( squares || iswithinHex(q,r,crH) )
				{
					if(th==0) //simple smoothing
					{
						{	pixels[ind] = (short)(filter.get_fval(q,r,r2) );} //replacement	
					}
					else
					{
						final int val = (int)pixels[ind];
						final int rep = filter.get_fval(q,r,r1); //detection
						if( filter_bright )
						{
							if( val - rep >= th)
							{	pixels[ind] = (short)(reuse ? rep : filter.get_fval(q,r,r2) );} //replacement
						}
						if( filter_dark )
						{
							if(rep - val >= th)
							{	pixels[ind] = (short)(reuse ? rep : filter.get_fval(q,r,r2) );} //replacement
						}
					}
					
				}
			}
		}
		return;
	}
	
	private boolean iswithinHex(final int q, final int r, final int cradius)
	{
		//if(q < 0 || r < 0)
		//{	IJ.error(getClass().getSimpleName() + ".iswithinHex: image coordinates cannot be negative");}
		//cube coords of target hex
		final int px = q - (r - (r&1)) / 2;
		final int pz = r;
		final int py = -px - pz;
		//offset hex in cube coords
		final int dx = px - cx;
		final int dy = py - cy;
		final int dz = pz - cz;
		return (dx >= -cradius && dx <  cradius &&
				dy >  -cradius && dy <= cradius &&
				dz >= -cradius && dz <  cradius );
	}
	
	private abstract class Filter
	{
		protected short[] old_pixels = null;
		public void setOldPix(final short[] pix)
		{
			old_pixels = Arrays.copyOf(pix,pix.length);
		}
		public abstract int get_fval(final int q, final int r, final int rH);
		protected int kernelsize = 0;
		protected int[] vals = null;
		protected int eval_fval() //no checking for sensible kernelsize, vals or mode
		{
			switch(mode)
			{
				case 0: //median
					Arrays.sort(vals);
					if(kernelsize%2 == 1)
					{	return vals[kernelsize/2];} //central value
					else if(kernelsize > 1)
					return (int)(0.5*(vals[kernelsize/2-1]+vals[kernelsize/2])+1);
				case 1: //mean
				{
					double sum = (vals[0]);
					for(int i = 1; i < kernelsize; ++i )
					{	sum += (vals[i]);}
					return (int)(sum/kernelsize + 0.5);
				}
				case 2: //min aka erode
				{
					int min = vals[0];
					for(int i = 1; i < kernelsize; ++i )
					{	
						if (vals[i] < min)
						{	min = vals[i];}	
					}
					return min;
				}
				case 3: //max aka dilate
				{
					int max = vals[0];
					for(int i = 1; i < kernelsize; ++i )
					{	
						if (vals[i] > max)
						{	max = vals[i];}	
					}
					return max;
				}
					
			}
			return(65535);
		} 
	}
	
	
	
	private class HexFilter extends Filter
	{
		public HexFilter() {};
		
		private boolean hexagonaltrap(final int[] cube, final int N)
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
				
		public int get_fval(final int q, final int r, final int rH)
		{
			final int px = q - r / 2;
			final int pz = r;
			final int py = -px - pz;
			kernelsize = weight;
			for(int i = 1; i <= rH; ++i)
			{	kernelsize += 6*i;} 
			vals = new int[kernelsize];
			int last_pix = 0; 
		
			for(int dx = -rH; dx <= rH; ++dx)
			{
				int dz1 = -rH;
				int dz2 = rH - dx;
				if ( dx < 0 )
				{
					dz1 = -rH - dx;
					dz2 = rH;
				}
				for(int dz = dz1; dz <= dz2; ++dz)
				{
					final int nx = px + dx ;
					final int nz = pz + dz ;
					final int ny = -nx - nz;
					
					int[] cube = {nx-cx, ny-cy, nz-cz };
					hexagonaltrap(cube, crH);
					
					final int nnx = cube[0] + cx;
					final int nnz = cube[2] + cz;
					
					final int nq = nnx + ( nnz - (nnz&1) ) / 2;
					final int nr = nnz;
					if( (dx != 0) || (dz != 0) )
					{	vals[last_pix++] = (old_pixels[nq + nr * imp_width] & 0xffff);}
					else for(int w = 0; w < weight; ++w)
					{	vals[last_pix++] = (old_pixels[nq + nr * imp_width] & 0xffff);}
				}
			}
			return eval_fval();
		}
	} //class HexFilter
	
	private class SquareFilter extends Filter
	{
		public SquareFilter() {};
		
		public int get_fval(final int q, final int r, final int rH)
		{
			kernelsize = (1+2*rH)*(1+2*rH) - 1 + weight;
			vals = new int[kernelsize];
			int last_pix = 0; 
			
			for(int y = r-rH; y <= r + rH; ++y)
			{
				int y0 = y % imp_height;
				if( y0 < 0)
				{	y0 += imp_height;}
				for(int x = q-rH; x <= q + rH; ++x)
				{
					int x0 = x % imp_width;
					if(x0 < 0) 
					{	x0 += imp_width;}
					if( (y != r) || (x != q) )
					{	vals[last_pix++] = (old_pixels[x0 + y0 * imp_width] & 0xffff);}
					else for (int w = 0; w < weight; ++w)
					{	vals[last_pix++] = (old_pixels[x0 + y0 * imp_width] & 0xffff);}	
				}
			}
			return eval_fval();
		}
	} //class SquareFilter
}
