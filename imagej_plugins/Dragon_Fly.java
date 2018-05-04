import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.FileSaver;
import java.awt.*;
import java.awt.geom.*;
import ij.plugin.*;
import java.util.Random;
import java.util.Scanner;
import java.util.TreeMap;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.IOException;

import java.nio.ByteBuffer;
import java.text.DecimalFormat;

public class Dragon_Fly implements PlugIn
{
	Random random = new Random();

	double hexedgelen = 6.0;
	int bondlength = 4;
	/* a1 vector is along zigzag direction */
	double a1x = hexedgelen * Math.sqrt(3);
	double a1y = 0.0;
	/* a2 vector is along armchair direction*/
	double a2x = 0.5 * a1x;
	double a2y = 1.5 * hexedgelen;
	/* the hexagonal lattice indices*/
	int na1 = 0;
	int na2 = 0;
	/* marks the middle of a square image*/
	//coords of active Hex
	int cq = 0;
	int cr = 0;
	int cx = cq - cr / 2;
	int cz = cr;
	int cy = -cx - cz;
	int cradius = 0;//only the hex at the center
	//cube coords of active translation
	int tx = 0;
	int tz = 0;
	int ty = -tx - tz;
	/* hexagonal lattice coordinates*/
	double na1f = 0.0;
	double na2f = 0.0;

	/* rectangular supercell*/
	double ax = hexedgelen * Math.sqrt(3);
	double axh = 0.5 * ax;
	double ay = 3.0 * hexedgelen;

	/* overlapping area in square pixels*/
	double overlap_area = 1.0;

	/* b1 and b2 are the regular "screen" basis bectors*/
	double b1x = 1.0;
	double b1y = 0.0;
	double b2x = 0.0;
	double b2y = 1.0;

	/*offset for sampling coordinates*/
	int sampleoffsetX = 0;
	int sampleoffsetY = 0;

	/*offsets for sampling the source image*/
	double pixeloffsetX = 0.0;
	double pixeloffsetY = 0.0;
	double alpha = 0.0;

	double old_pixeloffsetX = 0.0;
	double old_pixeloffsetY = 0.0;
	double old_alpha = 0.0;

	boolean setup_canceled = false;
	boolean limit_to_hex = false;
	boolean pseudo_hexes = false;
	boolean suppress_morie = false;
	boolean compensate_subpixels = true;
	boolean and_quit = false;
	boolean contrast_patternig = false;
	boolean full_hexagons = false;
	boolean split_image = false;
	boolean periodic_image = false;
	boolean crop_upon_rotation = false;
	boolean mark_center = false;
	boolean update = false;
	boolean append = false;
	boolean no_output = false;
	boolean suggest_newOffsets = false;
	boolean sample_and_view = false;
	boolean interactive_mode = false;
	boolean manual_scan_area = false;
	boolean use_mask = false;
	boolean resample_image = false;
	boolean random_sampling = true;
	boolean count_sampling = false;
	boolean use_old_bg = true;
	boolean lattice_sites = true;

	String splitname = "sf.tiff";
	String hexviewname = "hexview.tiff";
	String maskname = "<none>";
	String resamplename ="resampled.tiff";
	ImagePlus maskImg = null;
	
	int last_sfq = 0;
	int last_sfr = 0;
	
	int rotsteps = 1;
	int x_trans_inp = 0;
	int z_trans_inp = 0;
	//int Hexq,Hexr,Hexradius;

	int subframelength = 6;
	int viewing_border = 0;

	int sourcewidth = -1;
	int sourceheight = -1;
	int stacksize = 1;

	int outputwidth = -1;
	int outputheight = -1;

	int firstframe = 1;
	int lastframe = -1;
	int lastsubframe = -1;

	int hexcount = -1;
	int pixelcount = -1;

	int sfxstart = 0;
	int sfzstart = 0;
	int subpixels = 5;

	//variables for determining the center of hexareas
	long massX = 0;
	long massZ = 0;
	long massT = 0;
	double centX = 0;
	double centZ = 0;
	double suggestX = 0;
	double suggestY = 0;
	
	/*Range for subframe scanning*/
	// x and z of cube coords of corner points, must be convex polygon
	int point1_x = 0;
	int point1_z = 0;
	int point2_x = 0;
	int point2_z = 2;
	int point3_x = 2;
	int point3_z = 2;
	int point4_x = 2;
	int point4_z = 0;
	//cube normal vectors to the lines connecting points ij
	int norm12x = 0;
	int norm12z = 0;
	int norm23x = 0;
	int norm23z = 0;
	int norm34x = 0;
	int norm34z = 0;
	int norm41x = 0;
	int norm41z = 0;
	//constants for the lines along the edges
	int edge_D12 = 0;
	int edge_D23 = 0;
	int edge_D34 = 0;
	int edge_D41 = 0;
	
	int skprod(int v1x, int v1z, int v2x, int v2z)
	{
		int v1y = -v1x - v1z;
		int v2y = -v2x - v2z;
		return ( v1x*v2x + v1y*v2y + v1z*v2z);
	}
	
	//updates border normal vectores and the middle point from point1,2,3,4_x,z coords
	private void update_scan_frame()
	{
		int pointM_x = (point1_x + point2_x + point3_x + point4_x) / 4;
		int pointM_z = (point1_z + point2_z + point3_z + point4_z) / 4;
		//int pointM_y = -pointM_x - pointM_z;
		//edge from point1 to point2
		int v12x = point2_x - point1_x;
		int v12z = point2_z - point1_z;
		//edge from point2 to point3
		int v23x = point3_x - point2_x;
		int v23z = point3_z - point2_z;
		//edge from point3 to point4
		int v34x = point4_x - point3_x;
		int v34z = point4_z - point3_z;
		//edge from point4 to point1
		int v41x = point1_x - point4_x;
		int v41z = point1_z - point4_z;
		
		//simply assume one kind of normal vector for each edge
		norm12x =  2 * v12z + v12x;
		norm12z = -2 * v12x - v12z; 
		edge_D12 = skprod(norm12x,norm12z,point1_x,point1_z);
		/*
		if(edge_D12 != skprod(norm12x,norm12z,point2_x,point2_z))
		{
			IJ.error("error in edge_D12");
		}
		*/
		norm23x =  2 * v23z + v23x;
		norm23z = -2 * v23x - v23z; 
		edge_D23 = skprod(norm23x,norm23z,point2_x,point2_z);
		/*
		if(edge_D23 != skprod(norm23x,norm23z,point3_x,point3_z))
		{
			IJ.error("error in edge_D23");
		}
		*/
		norm34x =  2 * v34z + v34x;
		norm34z = -2 * v34x - v34z; 
		edge_D34 = skprod(norm34x,norm34z,point3_x,point3_z);
		/*
		if(edge_D34 != skprod(norm34x,norm34z,point4_x,point4_z))
		{
			IJ.error("error in edge_D34");
		}
		*/
		norm41x =  2 * v41z + v41x;
		norm41z = -2 * v41x - v41z; 
		edge_D41 = skprod(norm41x,norm41z,point4_x,point4_z);
		/*
		if(edge_D41 != skprod(norm41x,norm41z,point1_x,point1_z) )
		{
			IJ.error("error in edge_D41");
		}
		*/	
		//test if we have to negate the normal vectors	
		if(	skprod(norm12x,norm12z,pointM_x,pointM_z) < edge_D12 )
		{
			norm12x = -norm12x;
			norm12z = -norm12z;
			edge_D12 = skprod(norm12x,norm12z,point1_x,point1_z);
			/*
			if(	skprod(norm12x,norm12z,pointM_x,pointM_z) < edge_D12  )
			{
				IJ.error("error in flipping norm12");
			}
			*/ 
		}
		
		if(skprod(norm23x,norm23z,pointM_x,pointM_z) < edge_D23 )
		{
			norm23x = -norm23x;
			norm23z = -norm23z;
			edge_D23 = skprod(norm23x,norm23z,point2_x,point2_z);
			/*
			if(	skprod(norm23x,norm23z,pointM_x,pointM_z)  < edge_D23 )
			{
				IJ.error("error in flipping norm23");
			}
			*/ 
		}
		
		if(	skprod(norm34x,norm34z,pointM_x,pointM_z ) < edge_D34 )
		{
			norm34x = -norm34x;
			norm34z = -norm34z;
			edge_D34 = skprod(norm34x,norm34z,point3_x,point3_z);
			/*
			if(	skprod(norm34x,norm34z,pointM_x,pointM_z ) < edge_D34 )
			{
				IJ.error("error in flipping norm34");
			}
			*/ 
		}
		
		if(	skprod(norm41x,norm41z,pointM_x,pointM_z ) < edge_D41 )
		{
			norm41x = -norm41x;
			norm41z = -norm41z;
			edge_D41 = skprod(norm41x,norm41z,point4_x,point4_z);
			/*
			if(	skprod(norm41x,norm41z,pointM_x,pointM_z ) < edge_D41 )
			{
				IJ.error("error in flipping norm41");
			}
			*/ 
		}
		/*
		if( skprod(norm12x,norm12z,norm34x,norm34z) > 0 )
		{
				IJ.error("normalvectors of edge12 and edge34 are not pointing to one another");
		}
		
		if( skprod(norm23x,norm23z,norm41x,norm41z) > 0  )
		{
				IJ.error("normalvectors of edge23 and edge34 are not pointing to one another");
		}
		*/
		IJ.log("edge_D12: " + edge_D12 + "edge_D23: " + edge_D23 + "edge_D34: " + edge_D34 + "edge_D41: " + edge_D41);
	}
	
	private boolean point_inside_frame(int px, int pz)
	{
		return ( (skprod(norm12x,norm12z,px,pz) > edge_D12) &&
			     (skprod(norm23x,norm23z,px,pz) > edge_D23) &&
			     (skprod(norm34x,norm34z,px,pz) > edge_D34) &&
			     (skprod(norm41x,norm41z,px,pz) > edge_D41) );
	}
	
	
	private boolean inside_scan_frame(int px, int pz, int rH)
	{
		if(!manual_scan_area) return true;
		
		int lrH = rH;
		if(full_hexagons)
		{
			++lrH;
		}
		
		return point_inside_frame(px + lrH,pz) &&
			   point_inside_frame(px - rH,pz) &&
			   point_inside_frame(px,pz + lrH) &&
			   point_inside_frame(px,pz - rH) &&
			   point_inside_frame(px+lrH,pz-lrH) &&
	  		   point_inside_frame(px-rH,pz-rH);
	}
	
	private  boolean point_on_mask(int px, int pz)
	{
		int q = px + pz/2;
		int r = pz;
		int w = maskImg.getWidth();
		int h = maskImg.getHeight();
		if( (q < 0) || (q >= w) || (r < 0) || (r >= h) )
		{
			return false;
		}
		int val = maskImg.getPixel(q,r)[0];
		return( (val > 0) && (val <= 65535) );
	}
	
	
	private boolean inside_mask(int px,int pz,int rH) 
	{
		if(!use_mask) return true;
		
		int lrH = rH;
		if(full_hexagons)
		{
			++lrH;
		}
		
		return point_on_mask(px + lrH,pz) &&
			   point_on_mask(px - rH,pz) &&
			   point_on_mask(px,pz + lrH) &&
			   point_on_mask(px,pz - rH) &&
			   point_on_mask(px+lrH,pz-lrH) &&
	  		   point_on_mask(px-rH,pz-rH);
	}
	

	private void setImgSampling(double newpixeloffsetX, double newpixeloffsetY, double newalpha)
	{
		pixeloffsetX = newpixeloffsetX;
		old_pixeloffsetX = newpixeloffsetX;
		pixeloffsetY = newpixeloffsetY;
		old_pixeloffsetY = newpixeloffsetY;
		alpha = newalpha;
		old_alpha = newalpha;
		b1x = Math.cos(alpha);
		b1y = -Math.sin(alpha);
		b2x = Math.sin(alpha);
		b2y = Math.cos(alpha);
	}
	
	private void setImgViewing()
	{
		pixeloffsetX = 0.0;
		pixeloffsetY = 0.0;
		alpha = 0.0;
		b1x = 1.0;
		b1y = 0.0;
		b2x = 0.0;
		b2y = 1.0;
	}
	
	private void restoreImgSampling()
	{
		setImgSampling(old_pixeloffsetX, old_pixeloffsetY, old_alpha);
	}
	

	//squares to "odd-r" hexes. The origin is offset
	// see http://www.redblobgames.com/grids/hexagons/
	private void putSquare(double nnb1, double nnb2)
	{
		na1f = (nnb1 * b1x + nnb2 * b2x + pixeloffsetX) / ax;
		na2f = (nnb2 * b2y + nnb1 * b1y + pixeloffsetY) / ay;

		//algorithm using a face centered rectangular unitcell
		/* Use the Voroni property, check which corner or if
		 * the center is closest
		 */

		na1 = (int) Math.floor(na1f);
		na2 = (int) Math.floor(na2f);

		double na1d = (na1f - (double) na1) * ax;
		double na2d = (na2f - (double) na2) * ay;

		double d0 = (na1d-0.5*ax) * (na1d-0.5*ax) +  (na2d-0.5*ay) * (na2d-0.5*ay); //center
		double d1 = na1d * na1d + na2d * na2d; //top left
		double d2 = (na1d-ax) * (na1d-ax) +  (na2d) * (na2d); //top right
		double d3 = (na1d) * (na1d) +  (na2d-ay) * (na2d-ay); //lower left
		double d4 = (na1d-ax) * (na1d-ax) + (na2d-ay) * (na2d-ay); //lower right

		na2 *= 2;
		if      (d0 <= d1 && d0 <= d2 && d0 < d3 && d0 < d4) //center
		{
			//centered hexagons have odd rownumbers
			++na2;
		}
		else if (d1 <= d2 && d1 < d3 && d1 < d4) //top left
		{
			//noting to do
		}
		else if (d2 < d3 && d2 < d4) //top right
		{
			++na1;
		}
		else if (d3 < d4) //lower left
		{
			na2 +=2;
		}
		else //lower right
		{
			++na1;
			na2 +=2;
		}
	}

	//performs close to center test for "odd-r" hexes
	private boolean iswithinHex(int q, int r)
	{
		if(q < 0 || q < 0)
		{
			IJ.error("DragonFly.iswithinHex: image coordinates cannot be negative");
		}
		//cube coords of target hex
		int px = q - r / 2;
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

	public void sethexedgelen(double newhexedgelen)
	{
		if(newhexedgelen > 0)
		{
			hexedgelen = newhexedgelen;
			//axial basis
			a1x = hexedgelen * Math.sqrt(3);
			a1y = 0; //should be a const anyways
			a2x = 0.5 * a1x;
			a2y = 1.5 * hexedgelen;
			//rectangular basis
			ax = hexedgelen * Math.sqrt(3);
			axh = 0.5 * ax;
			ay = 3.0 * hexedgelen;
		}
		else
		{
			IJ.error("Dragonfly requires positive hexedgelen: "+newhexedgelen);
		}
	}

	public void setsquashedhexes(double extH, double extV)
	{
		if((extH <= 0) && (extV <= 0))
		{
			IJ.error("Dragonfly requires positive width,height: " + extH + "," + extV);
		}
		//axial coordinates
		a1x = extH;
		a2x = 0.5 * a1x;
		a2y = 0.75 * extV;
		//rectangular coordinates
		ax = extH;
		axh = 0.5 * ax;
		ay = 1.5 * extV;

	}

	public void setAutocenter(int a1max, int a2max)
	{
		if(a1max < 0 || a2max < 0)
		{	IJ.error("DragonFly.setAutocenter: image limits cannot be negative");}
		if(full_hexagons)
		{
			--a1max;
			--a2max;
			if(a1max < 0 || a2max < 0)
			{
				IJ.error("DragonFly.setAutocenter: image is too small for a full hexagon");
				a1max = 1;
				a2max = 1;
			}
		}
		int centera1 = (a1max)/2;
		int centera2 = (a2max)/2;
		int center = centera2;
		if(centera1 < centera2)
		{	center = centera1;}
		//cube coords of the center
		cq = center;
		cr = center;
		cx = center - center / 2;
		cz = center;
		cy = -cx - cz;
		cradius = center;
	}

	public void setHexcenter(int q, int r, int rH)
	{
		//in principle negative centers are allowed
		cq = q;
		cr = r;
		cradius = rH;
		cx = q - r / 2;
		cz = r;
		cy = -cx - cz;
	}

	//Warning here we use axial coordinates
	public void setTranslation(int q, int r)
	{
		tx = q;
		tz = r;
		ty = -tx - tz;
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

	public void imgSampling( short[] sourcepix, short[] outpix)
	{
		double[] dbloutpix =  new double[outpix.length];;
		short[] samplecount = new short[outpix.length];
		double subfactor = 1.0/subpixels;
		double pixarea =  (double)sourcepix.length * subpixels * subpixels / (double)outpix.length;
		if (pixarea == 0) {pixarea = 1;}
		int subarea = subpixels * subpixels;
		for(int i = 0; i < sourcepix.length; ++i)
		{
			int x = i % sourcewidth;
			int y = i / sourcewidth;
			
			for(int subx = 0; subx < subpixels; ++subx)
			{
				for(int suby = 0; suby < subpixels; ++suby)
				{
						int runs = count_sampling ? sourcepix[i] : 1;
						int cts = count_sampling ? 1 : sourcepix[i];
						for(int k = 0; k < runs; ++k)
						{
							//take either top left corner or a random position inside the subpixel
							double posX = x + subfactor * (subx + (random_sampling ? random.nextDouble() : 0) );
							double posY = y + subfactor * (suby + (random_sampling ? random.nextDouble() : 0) );
							putSquare(posX,posY); //sets na1 and na2 to hexagonal q,r
							int sX = na1 + sampleoffsetX;
							int sY = na2 + sampleoffsetY;
							int sInd = sX + sY * outputwidth;
							if( (!limit_to_hex) || iswithinHex(sX,sY) )
							{
								if(k==0)
								{	++samplecount[sInd];}
								dbloutpix[sInd] += cts;
							}
						}
				}
			}
				
		}
		if(suppress_morie)
		{
			for(int i = 0; i < dbloutpix.length; ++i)
			{
				if(samplecount[i] > 0)
				{
					dbloutpix[i] = Math.round( dbloutpix[i]*pixarea / ((double)samplecount[i] ));
				}
			}
		}
		sfxstart = 1; //set default for later subframe scanning
		sfzstart = 1;
		double scaling = compensate_subpixels ? (subpixels*subpixels) : 1.0;
		for(int i = 0; i < outpix.length; ++i)
		{
			outpix[i] = (short) Math.round( ((double)dbloutpix[i]/(scaling)));
		}	
	}

	public int sfCounting(short[] sourcepix)
	{
		int rH = subframelength;
		int lrH = rH;
		if(full_hexagons)
		{
			++lrH;
		}
		//new frame coords
		int xnf = rH - rH / 2;
		int znf = rH;

		int x0 = 0;
		int z0 = 1 + 2*rH+1;

		if(!full_hexagons)
		{
			//the "0th" row and col may contain cut intensities
			int qq = rH;
			int rr = rH;
			x0 = qq - rr / 2;
			z0 = rr;
		}

		x0 += sfxstart;
		z0 += sfzstart;

		int q0 = x0 + z0 / 2;
		int r0 = z0;
		int picked_hexes = 0;
		int first_sfx = 0;
		int sfz = 0;
		boolean sfx_picked = false;
		int	last_sfx = (sourcewidth + sourceheight)/(2*rH);
		int min_sfz = (sourceheight)/(2*rH);
		last_sfq = 0;
		last_sfr = 0;
		do //scanning sfz
		{
			sfx_picked = false;
			int sfx = first_sfx;
			do //scanning sfx
			{
				int x2 = x0;
				int z2 = z0;
				if(full_hexagons)
				{
					x2 += sfx*(2*rH + 1) - sfz*(rH);
					z2 += sfz*(2*rH + 1) - sfx*(rH + 1);
				}
				else
				{
					x2 +=  sfx * 2 * rH - sfz * rH;
					z2 +=  sfz * 2 * rH - sfx * rH;
				}
				//cube -> odd-r
				int q = (x2 + z2 / 2) ;
				int r = z2;

				if(q >= rH && q <= sourcewidth-lrH &&
				   r >= rH && r <= sourceheight-lrH && 
				   inside_scan_frame(x2,z2,rH) && 
				   (!use_mask || inside_mask(x2,z2,rH) ) ) //we are fully inside the source img
				{
					if(! sfx_picked)
					{	first_sfx = sfx - 1;} //start searching in next line one to the left
					sfx_picked = true;
					++picked_hexes; //counting per frame


					if( (lastsubframe >= 1) && (picked_hexes >= lastsubframe))
					{   //the first test is redundant
						return picked_hexes;
					}
					int sfq = q / (2*rH);
					int sfr = r / (2*rH);
					last_sfq = (sfq > last_sfq) ? sfq : last_sfq;
					last_sfr = (sfr > last_sfr) ? sfr : last_sfr;
				}
				else
				{
					if(sfx_picked || sfx > last_sfx) // stop searching sfx
					{
						break;
					}
				}
				++sfx;
			} while(true);//loop sfx
			++sfz;
		} while (sfx_picked || sfz <= min_sfz); //loop sfz while there are any useful sfx found
		++last_sfq;
		++last_sfr;
		return picked_hexes;
	}

	public void imgScanning(short[] sourcepix, short[] outpix, short[] resamplepix, ImagePlus simg)
	{
		int rH = subframelength;
		int lrH = rH;
		int dim = 2*rH;
		if(full_hexagons)
		{
			++dim;
			++lrH;
		}
		//new frame coords
		int xnf = rH - rH / 2;
		int znf = rH;

		int x0 = 0;
		int z0 = 1 + 2*rH+1;

		if(!full_hexagons)
		{
			//the "0th" row and col may contain cut intensities
			int qq = rH;
			int rr = rH;
			x0 = qq - rr / 2;
			z0 = rr;
		}

		x0 += sfxstart;
		z0 += sfzstart;

		int q0 = x0 + z0 / 2;
		int r0 = z0;
		int picked_hexes = 0;
		int first_sfx = 0;
		int sfz = 0;
		boolean sfx_picked = false;
		int	last_sfx = (sourcewidth + sourceheight)/(2*rH);
		int min_sfz = (sourceheight)/(2*rH);
		do //scanning sfz
		{
			sfx_picked = false;
			int sfx = first_sfx;
			do //scanning sfx
			{
				int x2 = x0;
				int z2 = z0;
				if(full_hexagons)
				{
					x2 += sfx*(2*rH + 1) - sfz*(rH);
					z2 += sfz*(2*rH + 1) - sfx*(rH + 1);
				}
				else
				{
					x2 +=  sfx * 2 * rH - sfz * rH;
					z2 +=  sfz * 2 * rH - sfx * rH;
				}
				//cube -> odd-r
				int q = (x2 + z2 / 2);
				int r = z2;

				if(q >= rH && q <= sourcewidth-lrH &&
				   r >= rH && r <= sourceheight-lrH &&
				   inside_scan_frame(x2,z2,rH) && 
				   (!use_mask || inside_mask(x2,z2,rH) )) //we are fully inside the source img
				{
					double f = 0.7 + 0.6*random.nextDouble();
					if(! sfx_picked)
					{	first_sfx = sfx - 1;} //start searching in next line one to the left
					sfx_picked = true;
					++hexcount; // total counting
					++picked_hexes; //counting per frame
					int sfq = q/(2*rH);
					int sfr = r/(2*rH);
					short[] pixels = null;
					if(simg != null)
					{
						ImageStack simgSt = simg.getStack();
						pixels = (short[])simgSt.getPixels(hexcount);
					}
					short wsum = 0;
					for(int dx = -rH; dx <= rH ; ++dx)
					{
						for(int dz = Math.max(-rH,-rH-dx); dz <= Math.min(rH, rH-dx); ++dz)
						{
							int dy = -dx - dz;
							if(!full_hexagons && (dz==rH || dy==-rH || dx==rH) )
							{
								continue;
							}

							int px2 = x2 + dx;
							int pz2 = z2 + dz;

							//cube -> odd-r
							int q2 = (px2 + pz2 / 2);
							int r2 = pz2;
							short w = sourcepix[q2 + r2 * sourcewidth];
							wsum += w;
							if (pixels != null)
							{
								//newframe coords
								int pxnf = xnf + dx;
								int pznf = znf + dz;
								int qnf = (pxnf + pznf / 2);
								int rnf = pznf;
								if((qnf >= 0) && (rnf >= 0) && (qnf < dim) && (rnf < dim) )
								{
									pixels[qnf + rnf * dim] = w;
								}
								else
								{
									System.out.println("DragonFly.ImgScanning Range ERROR q,r: " + qnf + "," + rnf);
								}
							}
							if(!no_output)
							{
								boolean is_center = (dx == 0 && dz == 0 );
								if(!is_center || !mark_center )
								{
									if(contrast_patternig)
									{
										outpix[q2 + r2 * outputwidth] = (short) (1000*f);//(f*sourcepix[q2 + r2 * sourcewidth]);
									}
									else
									{
										outpix[q2 + r2 * outputwidth] = w;
									}
									++pixelcount;
								}
							}
						}
					}
					if(resample_image)
					{
						resamplepix[sfq + sfr * last_sfq] = wsum;
					}
					
					if( (lastsubframe >= 1) && (picked_hexes >= lastsubframe))
					{   //the first test is redundant
						return;
					}
				}
				else
				{
					if(sfx_picked || sfx > last_sfx) // stop searching sfx
					{
						break;
					}
				}
				++sfx;
			} while(true);//loop sfx
			++sfz;
		} while (sfx_picked || sfz <= min_sfz); //loop sfz while there are any useful sfx found
	}

	public void imgViewing(short[] sourcepix, short[] outpix)
	{
		for(int i = 0; i < outpix.length; ++i)
		{
			int sx = i % outputwidth;
			int sy = i / outputwidth;

			putSquare( sx, sy);
			int q = na1;
			int r = na2;
			if(viewing_border > 0)
			{
				//odd-r -> cube
				int cubeX = q - r / 2;
				int cubeZ = r;
				cubeX -= viewing_border;
				cubeZ -= viewing_border;
				//cube -> odd-r
				q = cubeX + cubeZ / 2;
				r = cubeZ;
			}
			else if (periodic_image)
			{
				//no border = periodic boundary conditions
				q = q % sourcewidth;
				r = r % sourceheight;
				if(q < 0) {q += sourcewidth;}
				if(r < 0) {r += sourceheight;}
			}

			if(q >= 0 && q < sourcewidth &&
				   r >= 0 && r < sourceheight)
			//hexpixel lies inside source img
			{
				if( (!limit_to_hex) || iswithinHex(q,r) )
				{
					outpix[i] += sourcepix[ q + r * sourcewidth];
				}
			}
		}
	}

	public void imgRotating(short[] sourcepix, short[] outpix)
	{
		if(use_old_bg)
		{	setImgBg(sourcepix, outpix);}
		for(int i = 0; i < sourcepix.length; ++i)
		{
			int q = i % sourcewidth;
			int r = i / sourcewidth;
			if( iswithinHex(q,r))
			{
				//cube coords of target hex
				int px = q - r / 2;
				int pz = r;
				int py = -px - pz;

				//offset hex in cube coords
				int dx = px - cx;
				int dy = py - cy;
				int dz = pz - cz;

				int[] dvec = {dx,dy,dz};
				//alternating sign
				int pm = -1;
				if ( (rotsteps&1) == 0)
				{	pm = 1;}
				int j = (0 + rotsteps) % 3;
				if(j < 0) { j += 3;}
				int rx = pm * dvec[j];
				j = (1 + rotsteps) % 3;
				if(j < 0) { j += 3;}
				int ry = pm * dvec[j];
				j = (2 + rotsteps) % 3;
				if(j < 0) { j +=3 ;}
				int rz = pm * dvec[j];

				//translate back
				int[] cube = {rx , ry , rz};
				//apply periodic boundary conditions
				if(!full_hexagons) //full hexagons have rotation symmetry anyways
				{
					hexagonaltrap(cube, cradius);
				}
				px = cube[0]+ cx;
				py = cube[1]+ cy;
				pz = cube[2]+ cz;

				// cube to odd-r coords
				q = (px + pz / 2) % outputwidth;
				r = pz % outputheight;

				if(q < 0) {q += outputwidth;}
				if(r < 0) {r += outputheight;}
				//System.out.println("rot" + x0 + "," + y0 + " -> " + xp + "," + yp);
				outpix[q + r * outputwidth] = sourcepix[i];
			}
		}
	}

	public void imgPickHexArea(short[] sourcepix, short[] outpix)
	{
		final int k = 1; // use k = 2*bondlength + 1 to visit only equivalent lattice sites
		final int rH = cradius;
		for(int dx = -rH; dx <= rH ; dx+=k)
		{
			for(int dz = Math.max(-rH,-rH-dx); dz <= Math.min(rH, rH-dx); dz+=k)
			{
				int dy = -dx - dz;
				if(!full_hexagons && (dz==rH || dy==-rH || dx==rH) )
				{
					continue;
				}
				boolean is_center = (dx == 0 && dz == 0 );
				if(!is_center || !mark_center )
				{
					int px = dx + cx;
					//int py = dy + cy;
					int pz = dz + cz;
					int q = (px + pz / 2) % outputwidth;
					int r = pz % outputheight;
					//Maybe no periodic boundary conditions would be better behaviour?
					if(q < 0) {q += outputwidth;}
					if(r < 0) {r += outputheight;}
					outpix[q + r * outputwidth] = sourcepix[q + r * sourcewidth];
				}
			}
		}
	}

	public void imgTranslating(short[] sourcepix, short[] outpix)
	{
		if(use_old_bg)
		{	setImgBg(sourcepix, outpix);}
		final int rH = cradius;
		for(int dx = -rH; dx <= rH ; dx++)
		{
			for(int dz = Math.max(-rH,-rH-dx); dz <= Math.min(rH, rH-dx); ++dz)
			{
				int dy = -dx - dz;
				if(!full_hexagons && (dz==rH || dy==-rH || dx==rH) )
				{
					continue;
				}

				int px = dx + cx;
				//int py = dy + cy;
				int pz = dz + cz;

				//Maybe no periodic boundary conditions would be better behaviour?
				int q = (px + pz / 2) % sourcewidth;
				int r = pz % sourceheight;
				if(q < 0) {q += sourcewidth;}
				if(r < 0) {r += sourceheight;}
				short w = sourcepix[q + r * sourcewidth];

				int[] cube = {dx+tx, dy+ty, dz+tz };
				hexagonaltrap(cube, rH);

				int npx = cube[0] + cx;
				//int npy = cube[1] + cy;
				int npz = cube[2] + cz;
				if(suggest_newOffsets)
				{
					hexagonaltrap(cube, 3*bondlength);
					if(Math.abs(cube[0]) + Math.abs(cube[1]) + Math.abs(cube[2]) < 3*bondlength / 2)
					//full hexagons inside the truncated shape just not including the atom centers
					{
						//calculate center of mass
						massT += w;
						massX += (cube[0] * w);
						//massY += (cube[1] * w);
						massZ += (cube[2] * w);
					}
				}


				//Maybe no periodic boundary conditions would be better behaviour?
				int nq = (npx + npz / 2) % outputwidth;
				int nr = npz % outputheight;
				if(nq < 0) {nq += outputwidth;}
				if(nr < 0) {nr += outputheight;}
				outpix[nq + nr * outputwidth] = w;
			}
		}
	}

	public void imgMirroring(short[] sourcepix, short[] outpix)
	{
		if(use_old_bg)
		{	setImgBg(sourcepix, outpix);}
		final int rH = cradius;
		for(int dx = -rH; dx <= rH ; dx++)
		{
			for(int dz = Math.max(-rH,-rH-dx); dz <= Math.min(rH, rH-dx); ++dz)
			{
				int dy = -dx - dz;
				if(!full_hexagons && (dz==rH || dy==-rH || dx==rH) )
				{
					continue;
				}

				//Mirrored along x-axis
				int[] cube = { dx, dz, dy };

				hexagonaltrap(cube, rH);

				int px = dx + cx;
				int pz = dz + cz;

				int npx = cube[0] + cx;
				int npz = cube[2] + cz;

				//Maybe no periodic boundary conditions would be better behaviour?
				int q = (px + pz / 2) % sourcewidth;
				int r = pz % sourceheight;
				if(q < 0) {q += sourcewidth;}
				if(r < 0) {r += sourceheight;}

				//Maybe no periodic boundary conditions would be better behaviour?
				int nq = (npx + npz / 2) % outputwidth;
				int nr = npz % outputheight;
				if(nq < 0) {nq += outputwidth;}
				if(nr < 0) {nr += outputheight;}

				outpix[nq + nr * outputwidth] = sourcepix[q + r * sourcewidth];
			}
		}
	}

	public void setImgBg(short[] sourcepix, short[] outpix)
	{
		if(sourcepix.length != outpix.length)
		{
			IJ.error("Dimensions dont match, cannot set image background");
			return;
		}
		for(int i = 0; i < sourcepix.length; ++i )
		{
			outpix[i] = sourcepix[i];
		}
	}


	public ImagePlus imgTransform(ImagePlus sourceimg, String operation, String outname)
	{
		ImageProcessor sourceimgp = sourceimg.getProcessor();
		ImageStack sourceSt = sourceimg.getStack();

		if ( (sourceimg.getType() != ImagePlus.GRAY16))
		{
			IJ.error("DragonFly please provide a GRAY16 formated version");
			return sourceimg;
		}

		sourcewidth = sourceimg.getWidth();
		sourceheight = sourceimg.getHeight();
		int inputstacksize = sourceimg.getStackSize();

		if( (lastframe > inputstacksize) || (lastframe < 1) )
		{
			lastframe = inputstacksize;
		}
		if(firstframe < 1)
		{
			firstframe = 1;
		}


		stacksize = 1 + lastframe - firstframe;

		int outputframes = 1 + lastframe - firstframe;
		boolean sampling = operation.equals("sample Hexagons");
		boolean scanning = operation.equals("scan Subframes");
		boolean viewing = operation.equals("view Hexagons");
		boolean pickHexArea = operation.equals("clone");
		boolean rotating = operation.equals("rotate");
		boolean translating = operation.equals("translate");
		boolean mirroring = operation.equals("mirror");
		//only one transformation may be done in a single call
		if(sampling)
		{
			//top left corner
			int a1min = 0;
			int a1max = 0;
			int a2min = 0;
			int a2max = 0;
			point1_x = 0;
			point1_z = 0;

			//top right corner
			if(subpixels == 1)
			{
				putSquare(sourcewidth, 0);
			}
			else
			{
				putSquare(sourcewidth + 1, 0);
			}
			a1min = (na1 < a1min) ? na1 : a1min;
			a1max = (na1 > a1max) ? na1 : a1max;
			a2min = (na2 < a2min) ? na2 : a2min;
			a2max = (na2 > a2max) ? na2 : a2max;
			point2_x = na1 - na2/2;
			point2_z = na2;

			//lower right corner
			if(subpixels == 1)
			{
				putSquare(sourcewidth, sourceheight);
			}
			else
			{
				putSquare(sourcewidth + 1, sourceheight + 1);
			}
			a1min = (na1 < a1min) ? na1 : a1min;
			a1max = (na1 > a1max) ? na1 : a1max;
			a2min = (na2 < a2min) ? na2 : a2min;
			a2max = (na2 > a2max) ? na2 : a2max;
			point3_x = na1 - na2/2;
			point3_z = na2;

			//lower left corner
			if(subpixels == 1)
			{
				putSquare( 0 , sourceheight);
			}
			else
			{
				putSquare( 0 , sourceheight + 1);
			}
			a1min = (na1 < a1min) ? na1 : a1min;
			a1max = (na1 > a1max) ? na1 : a1max;
			a2min = (na2 < a2min) ? na2 : a2min;
			a2max = (na2 > a2max) ? na2 : a2max;
			point4_x = na1 - na2/2;
			point4_z = na2;

			outputwidth  = a1max - a1min + 1;
			outputheight = a2max - a2min + 1;
			sampleoffsetX = -a1min;
			sampleoffsetY = -a2min;
			int dx = -a1min + a2min/2;
			int dz = -a2min;
			point1_x += dx;
			point2_x += dx;
			point3_x += dx;
			point4_x += dx;
			point1_z += dz;
			point2_z += dz;
			point3_z += dz;
			point4_z += dz;
			update_scan_frame();

			setAutocenter(outputwidth, outputheight);

		}
		else if(viewing)
		{
			outputwidth = (int)Math.ceil( (sourcewidth + 2 * viewing_border) * a1x );
			outputheight = (int)Math.ceil( (sourceheight + 2 * viewing_border) * a2y );
			setAutocenter(sourcewidth, sourceheight);
			setImgViewing();
		}
		else if(rotating || pickHexArea || translating || scanning || mirroring)
		{
			outputwidth = sourcewidth;
			outputheight = sourceheight;
			if(scanning && use_mask)
			{
				maskImg = WindowManager.getImage(maskname);
				if(maskImg == null)
				{
					IJ.error("Could not find " + maskname);
					use_mask = false;
				}	
			}
		}
		else
		{
			IJ.error("DragonFly.imgTransform: unrecognized image transformation: " + operation);
			return sourceimg;
		}
		ImagePlus output = null;
		ImageProcessor outputimgp = null;
		ImageStack outputSt = null;
		if(!update && !append && !no_output)
		{
			output = NewImage.createShortImage(outname, outputwidth , outputheight , stacksize , NewImage.FILL_BLACK);
			outputimgp = output.getProcessor();
			outputSt = output.getStack();
		}

		ImagePlus simg = null;
		ImagePlus resampleImg = null;
		ImageStack resampleSt = null;
		if(scanning)
		{
			short[] pixels = (short[]) sourceSt.getPixels(firstframe);
			int sfCount = outputframes * sfCounting(pixels);
			if(split_image)
			{
				int dim = 2*subframelength;
				if(full_hexagons) {++dim;}
				simg = NewImage.createShortImage(splitname, dim, dim , sfCount, NewImage.FILL_BLACK);
				simg.show();
				IJ.log("slicing into subframes " + splitname);
			}
			if(resample_image)
			{
				resampleImg = NewImage.createShortImage(resamplename, last_sfq, last_sfr, outputframes, NewImage.FILL_BLACK);
				resampleSt = resampleImg.getStack();
				resampleImg.show();
				IJ.log("resampling to " + resamplename);
			}
		}
		hexcount = 0;
		pixelcount = 0;
		massX = 0;
		massZ = 0;
		massT = 0;

		try
		{
			int so = 0;
			for(int ss = firstframe; ss <= lastframe; ++ss )
			{
				IJ.showProgress(ss-firstframe , lastframe - firstframe + 1 );
				++so;
				short[] sourcepix = null;
				short[] resamplepix = null;
				short[] outpix = null;
				
				
				if(update)
				{
					outpix = (short[]) sourceSt.getPixels(ss);
					sourcepix = new short[outpix.length];
					for(int i = 0; i < outpix.length; ++i )
					{
						sourcepix[i] = outpix[i];
					}
				}
				else if(append)
				{

					sourceimg.setSlice( sourceimg.getStackSize() );
					IJ.run(sourceimg,"Add Slice","");
					sourceSt = sourceimg.getStack();
					sourcepix = (short[]) sourceSt.getPixels(ss);
					outpix = (short[]) sourceSt.getPixels(sourceimg.getStackSize());
				}
				else if(no_output)
				{
					sourcepix = (short[]) sourceSt.getPixels(ss);
				}
				else
				{
					sourcepix = (short[]) sourceSt.getPixels(ss);
					outpix = (short[]) outputSt.getPixels(so);
				}
				if(sampling)
				{
					imgSampling(sourcepix, outpix);
				}
				else if(scanning)
				{
					if(resample_image)
					{
						resamplepix = (short[]) resampleSt.getPixels(ss);
					}
					imgScanning(sourcepix, outpix, resamplepix, simg);
					if( (ss == stacksize) && (!no_output) )
					{
						int total_pixels = sourcewidth * sourceheight * stacksize;
						IJ.log("" + hexcount + " hexagons a " + ((double)pixelcount/hexcount) + "  coverage " + ((double)pixelcount/total_pixels) );
					}
				}
				else if(viewing)
				{
					imgViewing(sourcepix, outpix);
				}
				else if (rotating)
				{
					imgRotating(sourcepix, outpix);
				}
				else if(pickHexArea)
				{
					imgPickHexArea(sourcepix, outpix);
				}
				else if(translating)
				{
					imgTranslating(sourcepix, outpix);
				}
				else if(mirroring)
				{
					imgMirroring(sourcepix, outpix);
				}
			}
		}
		catch (Exception e)
        {
            e.printStackTrace();
            IJ.error("exception in DragonFly.imgTransform");
        }
        if(viewing) //we could simply always do that by default
        {	restoreImgSampling();}
        
        if(suggest_newOffsets && massT != 0)
        {
			suggest_newOffsets = false; //we dont make another popup default
			centX = (double)massX / massT;
			//centY = (double)massY / massT;
			centZ = (double)massZ / massT;
			suggestX = centX * a1x + centZ * a2x; //Maybe - here?
			suggestY = centX * a1y + centZ * a2y;
			IJ.log("Center of mass (cube) x,z " + centX + "," + centZ);
			IJ.log("Center of mass (pix) x,y " + suggestX + "," + suggestY);
			suggestX = old_pixeloffsetX - suggestX;
			suggestY = old_pixeloffsetY - suggestY;
			NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Please Confirm the new offsets");
			gd.addMessage("Set the new default values for hex sampling");
			gd.addMessage("previous offsets\nhorizontal " + old_pixeloffsetX + "\nvertical " + old_pixeloffsetY);

			gd.addNumericField("tilt", old_alpha, 6);
			gd.addNumericField("horizontal offset:", suggestX, 6);
			gd.addNumericField("vertical offset:", suggestY, 6);
			gd.showDialog();
			if( !gd.wasCanceled() )
			{
				double newalpha = gd.getNextNumber();
				double newpixeloffsetX = gd.getNextNumber();
				double newpixeloffsetY = gd.getNextNumber();
				setImgSampling(newpixeloffsetX, newpixeloffsetY, newalpha);
				IJ.log("updated sampling offset: " + old_pixeloffsetX + "," + old_pixeloffsetY + ",  alpha: " + old_alpha);
			}

        }
        if(!update && !append && !no_output)
        {
			return output;
		}
		else
		{
			no_output = false; //reset to default after operation
			return sourceimg;
		}
	}

	public void setuphex(String operation)
	{
		double newhexedgelen = hexedgelen;
		boolean sampling = operation.equals("sample Hexagons");
		boolean viewing = operation.equals("view Hexagons");
		
		hexviewname = WindowManager.makeUniqueName(hexviewname);
		
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Setup Hexagons");
		gd.addMessage("for " + operation);
		gd.addNumericField("hexedgelen:", newhexedgelen, 6,9,"pix");
		if(sampling)
		{
			gd.addNumericField("tilt:", old_alpha, 6,9,"pix");
			gd.addNumericField("horizontal offset:", old_pixeloffsetX, 6,9,"pix");
			gd.addNumericField("vertical offset:", old_pixeloffsetY, 6,9,"pix");
		}
		gd.addCheckbox("Use skewed hexagons", pseudo_hexes);
		if(sampling)
		{
			gd.addNumericField("subpixels", subpixels,0);
			gd.addCheckbox("Morie suppression", suppress_morie);
			gd.addCheckbox("Compensate subpixels", compensate_subpixels);
			gd.addCheckbox("Random sub position", random_sampling);
			gd.addCheckbox("individual counts (expensive)", count_sampling);
			gd.addCheckbox("view hexes (interactive)", sample_and_view);
			gd.addStringField("view: ", hexviewname, 30);
		}
		else if(viewing)
		{
			gd.addNumericField("border width", viewing_border,0);
			gd.addCheckbox("periodic image", periodic_image);
		}


		gd.showDialog();
		if (gd.wasCanceled())
		{
			setup_canceled = true;
			return;
		}
		newhexedgelen = gd.getNextNumber();
		if(sampling)
		{
			double newalpha = gd.getNextNumber();
			//newalpha = Math.abs(newalpha) % (Math.PI/6); //no need to restrain the angle
			double newpixeloffsetX = gd.getNextNumber();
			double newpixeloffsetY = gd.getNextNumber();
			setImgSampling(newpixeloffsetX, newpixeloffsetY, newalpha);
		}
		pseudo_hexes = gd.getNextBoolean();
		if(sampling)
		{
			subpixels = (int)gd.getNextNumber();
			if (subpixels < 1) {	subpixels = 1;}
			suppress_morie = gd.getNextBoolean();
			compensate_subpixels = gd.getNextBoolean();
			random_sampling = gd.getNextBoolean();
			count_sampling = gd.getNextBoolean();
			sample_and_view = gd.getNextBoolean();
			interactive_mode = sample_and_view;
			hexviewname = gd.getNextString();
			if(count_sampling && !random_sampling)
			{
				IJ.showMessage("Sampling counts individually, is only meaningfull with randomized sampling positions");
				count_sampling = false;
			}
		}
		else if(viewing)
		{
			viewing_border = (int)gd.getNextNumber();
			if(viewing_border < 0)
			{	viewing_border = 0;}
			periodic_image = gd.getNextBoolean();
		}
	
		if(compensate_subpixels || suppress_morie || random_sampling)
		{
			IJ.log("Morie suppression: " + suppress_morie + 
			", Compensate subpixels: " + compensate_subpixels +
			", Random sub position: " + random_sampling + 
			", individual counts: " + count_sampling);
		}
		
		
		
		if (pseudo_hexes)
		{
			double nextH = newhexedgelen * Math.sqrt(3.0);
			double nextV = 2 * newhexedgelen;

			NonBlockingGenericDialog gd2 = new NonBlockingGenericDialog("Squashed Hexagons");
			gd2.addMessage("for " + operation);
			gd2.addMessage("Choose horizontal(zigzag) and vertical(armchair) diameters");
			gd2.addNumericField("horizontal diameter:", nextH, 4);
			gd2.addNumericField("vertical diameter:", nextV, 4);

			gd2.showDialog();
			if (gd2.wasCanceled())
			{
				setup_canceled = true;
				return;
			}

			nextH = gd2.getNextNumber();
			nextV = gd2.getNextNumber();
			IJ.log("squashed hexagons w x h: " + nextH + " x " + nextV);
			IJ.log("skewness:" + ( 0.5 * Math.sqrt(3.0) * nextV / nextH ) );
			setsquashedhexes(nextH, nextV);
		}
		else
		{
			sethexedgelen(newhexedgelen);
			IJ.log("edgelen = " + hexedgelen );
			IJ.log("hexagons w x h: " + ax + " x " + (2 * hexedgelen));
		}
		if(viewing_border != 0)
		{	IJ.log("viewing border "+ viewing_border +" hexagons");}

		//sampling offset is periodic in rectangular supercell
		pixeloffsetX = pixeloffsetX % ax;
		if(pixeloffsetX < 0.0)
		{	pixeloffsetX += ax; }
		pixeloffsetY = pixeloffsetY % ay;
		if(pixeloffsetY < 0.0)
		{	pixeloffsetY += ay; }
		if( (pixeloffsetX != 0.0) || (pixeloffsetY != 0.0) || (alpha != 0.0))
		IJ.log("sampling offset: " + pixeloffsetX + "," + pixeloffsetY + "\ttilt: " + alpha);

	}

	public void setuprot(ImagePlus sourceimg)
	{
		boolean autocenter = true;
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Setup Rotation");
		gd.addNumericField("60° steps counterclockwise:", rotsteps, 0);
		gd.addCheckbox("source as background", use_old_bg);
		gd.addCheckbox("Biggest possible Hexagon", autocenter);
		gd.showDialog();
		if (gd.wasCanceled())
		{
			setup_canceled = true;
			return;
		}
		rotsteps = (int) gd.getNextNumber();
		use_old_bg = gd.getNextBoolean();
		
		autocenter = gd.getNextBoolean();

		if(autocenter)
		{
			ImageProcessor sourceimgp = sourceimg.getProcessor();
			setAutocenter( sourceimgp.getWidth(), sourceimgp.getHeight());
			IJ.log("hexagonal area r,q " + cq + "," +  cr + "   x,y,z: " + cx + "," + cy + "," + cz + "  radius = " + cradius);
		}
		else
		{
			setupHexArea();
		}
		IJ.log ("60° rotation steps" + rotsteps);
	}

	public void setuptrans(ImagePlus sourceimg)
	{
		boolean autocenter = true;
		
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Setup Translation");
		gd.addMessage("Translation with Periodic Boundary Conditions");
		//gd.addMessage("x up-right, z downwards");
		gd.addNumericField("Movement x:", x_trans_inp, 0);
		gd.addNumericField("Movement z:", z_trans_inp, 0);
		gd.addNumericField("bondlength:", bondlength, 0);
		gd.addCheckbox("walk by lattice sites( no = hexpixels )", lattice_sites);
		gd.addCheckbox("source as background", use_old_bg);
		gd.addCheckbox("full hexagons ", full_hexagons);
		gd.addCheckbox("Biggest possible Hexagon", autocenter);
		gd.addCheckbox("Suggest new offsets for sampling", suggest_newOffsets);
		gd.showDialog();
		if (gd.wasCanceled())
		{
			setup_canceled = true;
			return;
		}
		x_trans_inp = (int) gd.getNextNumber();
		z_trans_inp = (int) gd.getNextNumber();
		int xl = x_trans_inp;
		int zl = z_trans_inp;
		
		bondlength = (int) gd.getNextNumber();
		lattice_sites = gd.getNextBoolean();
		use_old_bg = gd.getNextBoolean();
		full_hexagons = gd.getNextBoolean();
		autocenter = gd.getNextBoolean();
		suggest_newOffsets = gd.getNextBoolean();
		if(lattice_sites)
		{
			xl = 2*x_trans_inp - z_trans_inp;
			zl = 2*z_trans_inp - x_trans_inp;
			xl *= bondlength;
			zl *= bondlength;	
		}
		
		if(autocenter)
		{
			ImageProcessor sourceimgp = sourceimg.getProcessor();
			setAutocenter( sourceimgp.getWidth(), sourceimgp.getHeight());
			IJ.log("hexagonal area r,q " + cq + "," +  cr + "   x,y,z: " + cx + "," + cy + "," + cz + "  radius = " + cradius);
		}
		else
		{
			setupHexArea();
		}
		setTranslation(xl, zl);
		IJ.log("translation x,y,z: " + tx + "," + ty + "," + tz);
	}



	public void setupHexArea()
	{

		NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Setup HexArea");
		gd.addNumericField("Position row:", cq, 0);
		gd.addNumericField("Position col:", cr, 0);
		gd.addNumericField("radius:", cradius, 0);
		gd.addCheckbox("source as background", use_old_bg);
		gd.addCheckbox("mark center with zero",mark_center);
		gd.addCheckbox("full hexagons ", full_hexagons);
		
		gd.showDialog();
		if (gd.wasCanceled())
		{
			setup_canceled = true;
			return;
		}
		int q = (int) gd.getNextNumber();
		int r = (int) gd.getNextNumber();
		int radius = (int) gd.getNextNumber();
		use_old_bg = gd.getNextBoolean();
		mark_center = gd.getNextBoolean();
		full_hexagons = gd.getNextBoolean();
		
		
		setHexcenter(q, r, radius);
		IJ.log("hexagonal area r,q " + cq + "," +  cr + "   x,y,z: " + cx + "," + cy + "," + cz + "  radius = " + cradius);
	}

	public void setupscan (ImagePlus sourceimg )
	{
		int[] idArray = WindowManager.getIDList();
		String[] titleArray = new String[idArray.length + 1]; // titles of opened images
		titleArray[0] = "<none>";
		for (int i = 0; i < idArray.length; ++i)
		{
			titleArray[i+1] = WindowManager.getImage(idArray[i]).getTitle();
		}
		if( WindowManager.getWindow(maskname) == null)
		{
			maskname = titleArray[0];
		}				
		resamplename = WindowManager.makeUniqueName(resamplename);
		splitname = WindowManager.makeUniqueName(splitname);
		
		NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Setup Subframes");
		gd.addMessage("The patterning must match the lattice !\nuse offset to skip partially sampled hexes");
		gd.addNumericField("subframelength (Pixels on an edge):", subframelength, 0);
		gd.addNumericField("offset x (right)", sfxstart, 0);
		gd.addNumericField("offset z (downright)", sfzstart, 0);
		gd.addNumericField("last subframe (<1 All)", lastsubframe, 0);
		gd.addCheckbox("split subframes", split_image);
		gd.addStringField("split output", splitname, 30);
		gd.addCheckbox("add contrast to hexagons", contrast_patternig);
		gd.addCheckbox("mark center with zero",mark_center);
		gd.addCheckbox("full hexagons ", full_hexagons);
		gd.addCheckbox("no output", no_output);
		gd.addCheckbox("custom scan area", manual_scan_area);
		gd.addChoice("mask for scanning (GRAY16)", titleArray, maskname );
		gd.addCheckbox("subframe Intensities (aka resampling)", resample_image);
		gd.addStringField("resample output", resamplename, 30);
		
		gd.showDialog();
		if (gd.wasCanceled())
		{
			setup_canceled = true;
			return;
		}
		subframelength = (int) gd.getNextNumber();
		sfxstart = (int) gd.getNextNumber();
		sfzstart = (int) gd.getNextNumber();
		lastsubframe = (int)gd.getNextNumber();
		split_image = gd.getNextBoolean();
		splitname = gd.getNextString();
		contrast_patternig = gd.getNextBoolean();
		mark_center = gd.getNextBoolean();
		full_hexagons = gd.getNextBoolean();
		no_output = gd.getNextBoolean();
		manual_scan_area = gd.getNextBoolean();
		maskname = gd.getNextChoice();
		use_mask = !maskname.equals(titleArray[0]);
		resample_image = gd.getNextBoolean();
		resamplename = gd.getNextString();

		resamplename = WindowManager.makeUniqueName(resamplename);
		splitname = WindowManager.makeUniqueName(splitname);
		if(full_hexagons && resample_image)
		{	
			IJ.error("Warning resampling is not possible with full hexagons");
			resample_image = false;
		}
	
		int imgwidth = sourceimg.getWidth();
		int imgheight = sourceimg.getHeight();

		int area =  1;
		if(full_hexagons)
		{
			area = 3 * subframelength * (subframelength + 1);
		}
		else
		{
			area = 3 * subframelength * subframelength;
		}
		IJ.log("subframeradius = " + subframelength + " Pixels =" + area);
		IJ.log("offset x,z " + sfxstart + "," + sfzstart);
		IJ.log("extra contrast " + contrast_patternig + "   full hexagons " + full_hexagons + "   last subframe " + lastsubframe);
		
		int p1_q = point1_x + point1_z / 2;
		int p2_q = point2_x + point2_z / 2;
		int p3_q = point3_x + point3_z / 2;
		int p4_q = point4_x + point4_z / 2;
		
		if(manual_scan_area)
		{
			
			
			NonBlockingGenericDialog gd2 = new NonBlockingGenericDialog("Setup Scanbox");
			gd2.addMessage("Enter the colum and row of the outmost possible");
			gd2.addMessage("subframe centers.");
			gd2.addMessage("p1->p2->p3->p4->p1 must be a convex shape.");
			gd2.addNumericField("point1 col", p1_q, 0);
			gd2.addNumericField("point1 row", point1_z, 0);
			gd2.addNumericField("point2 col", p2_q, 0);
			gd2.addNumericField("point2 row", point2_z, 0);
			gd2.addNumericField("point3 col", p3_q, 0);
			gd2.addNumericField("point3 row", point3_z, 0);
			gd2.addNumericField("point4 col", p4_q, 0);
			gd2.addNumericField("point4 row", point4_z, 0);
			
			gd2.showDialog();
			if (gd2.wasCanceled())
			{
				setup_canceled = true;
				return;
			}
			
			p1_q = (int) gd2.getNextNumber();
			point1_z = (int) gd2.getNextNumber();
			p2_q = (int) gd2.getNextNumber();
			point2_z = (int) gd2.getNextNumber();
			p3_q = (int) gd2.getNextNumber();
			point3_z = (int) gd2.getNextNumber();
			p4_q = (int) gd2.getNextNumber();
			point4_z = (int) gd2.getNextNumber();
			
			point1_x = p1_q - point1_z/2;
			point2_x = p2_q - point2_z/2;
			point3_x = p3_q - point3_z/2;
			point4_x = p4_q - point4_z/2;
				
		}
		update_scan_frame();
		IJ.log("scan limits (col,row) P1: " + p1_q + "," + point1_z + "  P2: " + 
		p2_q + "," + point2_z + "  P3: " + p3_q + "," + point3_z + "  P4: " +
		p4_q + "," + point4_z);
	}

	public void run(String arg)
	{
		try
		{
			String operation = "view Hexagons";
			String outname = "hexview.tiff";

			boolean keep_original = true;
			do
			{
				NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Dragon_Fly");
				int[] idArray = WindowManager.getIDList(); // list of all opened images (IDs)
				if (idArray == null)
				{
					//TODO call Graphene_Generator.run() directly here
					try
					{
						IJ.run("Graphene Generator");
						idArray = WindowManager.getIDList();
						operation = "view Hexagons";
						outname = "Graphene.tiff";
						periodic_image = true;
					}
					catch(RuntimeException re)
					{
						if(re.getMessage().equals("Macro canceled"))
						{	//User has canceled Graphene Generator
							return;
						}
						else
						{
							throw re;
						}
					}
					if(idArray == null) //just to be on the safe side
					{
						IJ.error("no input available quiting Dragon Fly");
						return;
					}
				}
				outname = WindowManager.makeUniqueName(outname);
				String[] titleArray = new String[idArray.length]; // titles of opened images
				for (int i = 0; i < idArray.length; ++i)
				{
					titleArray[i] = WindowManager.getImage(idArray[i]).getTitle();
				}
				String[] Operations = {"sample Hexagons","scan Subframes", "view Hexagons", "clone" ,"rotate","translate","mirror"};
				String candidate = WindowManager.getCurrentImage().getTitle();
				if (candidate == null || candidate.equals(""))
				{	candidate = titleArray[titleArray.length-1];}

				gd.addMessage("Choose an image and a transformation");
				gd.addChoice("Source image (GRAY16)", titleArray, candidate );
				gd.addNumericField("First frame to process", firstframe, 0);
				gd.addNumericField("Last frame to process (<1 End)", lastframe, 0);
				gd.addChoice("Transformation", Operations, operation );
				gd.addMessage("You may \"<update>\" or \"<append>\" the source ");
				gd.addStringField("Output", outname, 30);
				gd.addCheckbox("detect largest Hexagon", limit_to_hex);
				gd.addCheckbox("quit after Transformation", and_quit);
				gd.addCheckbox("keep original image", keep_original);

				gd.showDialog();
				if (gd.wasCanceled())
				{
					return;
				}

				String sourcetitle = gd.getNextChoice();
				firstframe = (int)gd.getNextNumber();
				lastframe = (int)gd.getNextNumber();
				operation = gd.getNextChoice();
				outname = gd.getNextString();
				limit_to_hex = gd.getNextBoolean();
				and_quit = gd.getNextBoolean();
				keep_original = gd.getNextBoolean();

				update = outname.equals("<update>");
				append = outname.equals("<append>");
				boolean sampling = operation.equals("sample Hexagons");
				boolean scanning = operation.equals("scan Subframes");
				boolean viewing = operation.equals("view Hexagons");
				boolean pickHexArea = operation.equals("clone");
				boolean rotating = operation.equals("rotate");
				boolean translating = operation.equals("translate");
				boolean mirroring = operation.equals("mirror");

				if( (update || append) && (viewing || sampling ) )
				{
					IJ.error("<update>/<append> are not possible for :" + operation);
					continue;
				}

				do //interactive mode shortcuts the main menu
				{
					outname = WindowManager.makeUniqueName(outname);
					
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
						IJ.log("source:" + sourcetitle + "   frames " + firstframe + " to End");
					}
					else
					{
						IJ.log("source:" + sourcetitle + "   frames " + firstframe  +" to " + lastframe);
					}

					IJ.log("transformation: " + operation);

					ImagePlus sourceimg = WindowManager.getImage(sourcetitle);
					if(sourceimg == null)
					{
						IJ.error("Could not find " + sourcetitle);
						continue;
					}
					int oldstacksize = sourceimg.getStackSize();
					sourceimg.show();
					setup_canceled = false;
					if(limit_to_hex)
					{

						setAutocenter( sourceimg.getWidth(), sourceimg.getHeight());
						if(! (translating || scanning || rotating || mirroring))
						{   //output only if there would not be another query for confirmation
							IJ.log("hexagonal center: " + cq + "," +  cr + "  radius = " + cradius);
						}
					}

					if (sampling || viewing ) //we need to define hex pixels
					{
						setuphex(operation);
					}
					if(scanning)
					{
						setupscan ( sourceimg );
					}

					if(rotating)
					{
						setuprot( sourceimg );
					}

					if(translating)
					{
						setuptrans( sourceimg );
					}

					if(pickHexArea || mirroring)
					{
						setupHexArea();
					}
					if (setup_canceled)
					{
						IJ.error(operation + " canceled");
						interactive_mode = false;
						continue;
					}

					ImagePlus tmp = imgTransform(sourceimg, operation, outname);

					if(!update && !append && tmp != null)
					{
						tmp.show();
						IJ.log("created " + outname);
					}
					else
					{
						tmp.updateAndRepaintWindow();
						if(update)
						{
							IJ.log("updated " + sourcetitle + " frames " + firstframe + " to " + lastframe);
						}
						else //append
						{
							int newframes = tmp.getStackSize() - oldstacksize;
							IJ.log("appended " + newframes + " frames to " + sourcetitle);
							firstframe = oldstacksize + 1;
							lastframe = -1;
						}
					}
					
					if(sample_and_view)
					{
						ImagePlus tmp2 = imgTransform(tmp, "view Hexagons", hexviewname);
						if( tmp2 != null)
						{
							tmp2.show();
							IJ.log("created " + hexviewname);
						}
					}
				
					if(!keep_original && !update)
					{
						sourceimg.changes = false;
						sourceimg.close();
					}
					else if (!keep_original && update)
					{
						keep_original = true;
						IJ.log("DragonFly refuses user request to update AND discard " + sourcetitle);
					}
				}
				while (interactive_mode);

			} while(!and_quit);

		}
		catch (Exception e)
        {
            e.printStackTrace();
        }
	}
}
