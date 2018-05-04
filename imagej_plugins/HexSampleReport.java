import ij.*;
import java.util.Scanner;
import java.util.HashMap;
import java.util.regex.Pattern;

public class HexSampleReport
{
	double elapsed = -1.0;
	double hexedgelen = 1.0;
	double tilt = 0.0;
	double ellipseA = 1.0;
	double ellipseB = 1.0;
	double ellipsePhi = 0.0;
	double excent = 0.0;
	double offset_x = 0.0;
	double offset_y = 0.0;
	double merit = 0.0;
	double old_merit = 0.0;
	double hex_merit = 0.0;
	double old_hex_merit = 0.0;
	double coverage = 0.0;
	double free_area = 0.0;
	double raw_avg = 0.0; //the average intensity inside complete unitcells
	double raw_std = 0.0; 
	double hex_shape = -1000.0;
	double moment2 = 0.0;
	double hex_avg = 0.0;
	double hex_std = 0.0; 
	double hex_low = 0.0;
	double hex_high = 0.0;
	double contrast = 0.0;
	double positionQ = 0.0;
	double mirrorQ = 0.0;
	double shapeQ = 0.0;
	double tuning = -1.0;
	int offset_hexq = 1;
	int offset_hexr = 1;
	int resets = 0;
	int subframes = 0;
	int frNr = 0;
	int mkNr = 0;
	int state = 0;   //set and interpreted by Hex_Sampler
	int peaks = -1; //number of first Bragg reflexes
	int worker = -1;
	int runs = 0;
	int stability = 2;
	int group = -1; //Not even a valid group before not even defaults were provided
	double stageX = Double.NaN;
	double stageY = Double.NaN;
	double sym = -1.0;
	int stageN = 0;
	int imgID = -1;
	int repID = 1;
	
	boolean use_mirror = false;
	String label = null;
	public static double def_excent = 1.0;
	public static double def_phi = 0.0;
	public static int current_group = -1;
	ImagePlus Idata = null;
	ImagePlus Mask = null;
	
	public static boolean label_ok = true;
	private boolean unread_label = true;
	
	
	final static public String headerline = "#label\tslice\tgraphene\tpeaks\ttuning\ttilt\tellA\tellB\tellPhi" + 
			"\thel\toffX\toffY\tNsf\toffQ\toffR\ttime\tmerit\tuc_avg\thex_avg\tcontrast"+
			"\tpositionQ\tmirrorQ\thex_low\thex_std\tuc_std\texcent\tshapeQ\thex_high\thex_merit" +
			"\thex_shape\tmoment2\tstate\tstageX\tstageY\tstageN\timgID\tsym";
			
	static String keyline = headerline.substring(1);
	//make this public maybe we want to setup some external aliases
	static public HashMap<String,String> aliases;
	static
	{
		aliases = new HashMap<String,String> ();
		aliases.put("filename","label");
		aliases.put("numpeak","peaks");
		aliases.put("ella","ellA");
		aliases.put("ellb","ellB");
		aliases.put("ellphi","ellPhi");
		aliases.put("Hex_Shape","hex_shape");
		aliases.put("contrast(std/mean)","contrast");
	} 
	
	static public void set_key_line(String new_keys)
	{
		keyline = new_keys;
		++current_group;
	}
	
	static private void get_an_alias(String plain_key)
	throws  UnsupportedOperationException, 
			ClassCastException,
			NullPointerException, 
			IllegalArgumentException
	{
		String default_key = "<none>";
		String new_key = IJ.getString("provide an alias for \"" + plain_key +"\"", default_key ).trim();
		if( new_key.equals("") )
		{	new_key = default_key;}
		aliases.put(plain_key,new_key);
		IJ.log(plain_key + " -> " + new_key);
	}
	
	public HexSampleReport(){}; //default contructor somehow explicitely needed
	
	public HexSampleReport(double hel, double radius, double n_tilt, int n_frNr) //with default values
	{
		hexedgelen = hel;
		ellipseA = radius * Math.sqrt(def_excent);
		ellipseB = radius / Math.sqrt(def_excent);
		ellipsePhi = def_phi;
		tilt = n_tilt;
		frNr = n_frNr;
		mkNr = n_frNr; //resonable estimate will be adjusted anyways 
		repID = n_frNr;
		group = 0; //the default group, if at least defaults were provided
	}
	
	public String toString()
	{
		update_redundant();
		String line =      label +
					"\t" + frNr +
					"\t" + free_area +
					"\t" + peaks +
					"\t" + tuning +
					"\t" + tilt +
					"\t" + ellipseA +
					"\t" + ellipseB +
					"\t" + ellipsePhi +	
					"\t" + hexedgelen +
					"\t" + offset_x +
					"\t" + offset_y +
					"\t" + subframes +
					"\t" + offset_hexq +
					"\t" + offset_hexr +
					"\t" + elapsed +
					"\t" + merit +
					"\t" + raw_avg +
					"\t" + hex_avg +
					"\t" + contrast +
					"\t" + positionQ +
					"\t" + mirrorQ +
					"\t" + hex_low +
					"\t" + hex_std +
					"\t" + raw_std +
					"\t" + excent +
					"\t" + shapeQ +
					"\t" + hex_high +
					"\t" + hex_merit +
					"\t" + hex_shape +
					"\t" + moment2 +
					"\t" + state +
					"\t" + stageX +
					"\t" + stageY +
					"\t" + stageN + 
					"\t" + imgID +
					"\t" + sym;
		return line;
	}
	public boolean parseString(String line, ImagePlus Idata, ImagePlus Mask)
	{	return parseString(line, Idata, null, null, false, Mask);}
	
	public boolean parseString(String line, ImagePlus Idata, ImageStack IdataStO, ImageStack MaskStO, boolean reorder_Idata, ImagePlus Mask)
	{
		this.Idata = Idata; //keep a local reference
		this.Mask = Mask;   //keep a local reference
		group = current_group;
		String key = null;
		String plain_key = null;
		try
		{
			ImageStack IdataSt = Idata.getStack();
			
			int impDepth = IdataSt.getSize();
			ImageStack MaskSt = null;
			int mskDepth = 0;
			if(Mask != null) //Hex_Collector does not need the Mask
			{
				MaskSt = Mask.getStack();
				mskDepth = MaskSt.getSize();
			}
			if(reorder_Idata && IdataStO == null)
			{
				IJ.log("Error IdataStO cannot be null if sorting is requested");
				reorder_Idata = false;
			}
			double ea   = -1.0;	
			double eb   = -1.0;
			double ephi	= Double.NaN;
			Scanner ls = new Scanner(line);
			Scanner ks = new Scanner(keyline);
			Pattern ptrn = Pattern.compile("[\\s,]+"); //support for csv
			ls.useDelimiter(ptrn);
			ks.useDelimiter(ptrn);
			boolean slice_found = false;
			boolean mask_found = false;
			boolean has_label = false;
			while( ls.hasNext() && ks.hasNext() )
			{
				
				plain_key = ks.next().trim();
				key = plain_key;
				
				//We could get funny redirections here by looping ;)
				if(aliases.containsKey(plain_key))
				{	key = aliases.get(plain_key);}
				
				if(key.equals("label"))
				{
					label = ls.next();
					if(label.equals("null")) //That is how missing labels get written to file
					{
						label = null; //the literal meaning
					}
					if(label != null)
					{
						int lpos = label.lastIndexOf('.');
						if(lpos > 0)
						{
							String fext = label.substring(lpos);
							if( fext.startsWith(".tif") || fext.startsWith(".dm3")) //drop .dm3 , .tif or .tiff extensions silently
							{	label = label.substring(0, lpos);}
						}
						//look for matching slicelabel in Idata
						int fr = 1;
						//Here we assume either all or none are labelled
						has_label= true;
						while ( fr <= impDepth &&
								(IdataSt.getSliceLabel(fr) != null) &&
								!(IdataSt.getSliceLabel(fr).contains(label)) )
								// !( IdataSt.getSliceLabel(fr).startsWith(label) ) && 
								// !( label.startsWith(IdataSt.getSliceLabel(fr))) )
						{	fr++;}		
						slice_found = ( (fr <= impDepth) && (IdataSt.getSliceLabel(fr) != null) );
						if( (fr <= impDepth) && (IdataSt.getSliceLabel(fr)==null)) //first test is redundant
						{
							IdataSt.setSliceLabel(label, fr);
							Idata.setStack(IdataSt); //update this one immediately since it may be also be reordered
						}
						
						if( (fr <= impDepth) && (IdataSt.getSliceLabel(fr) != null) ) //both test are redundant
						{
								
							if ( reorder_Idata )
							{
								//handles to identified frame
								short[] pixels = (short[])IdataSt.getPixels(fr);
								String frlbl = IdataSt.getSliceLabel(fr);
								//put identified frame (fr) to target frNr
								IdataStO.setPixels(pixels, frNr );
								IdataStO.setSliceLabel(frlbl, frNr );
							}
							else
							{	frNr = fr;}
							
							
						}
						else 
						{	
							
							if(fr >= impDepth) 
							{	
								IJ.log("WARNING: could not find frame " + label); 
								//reply = false;	
							}
						}
						///Now look for matching label inside Mask
						//Here we assume either all or none are labelled
						if(Mask != null) //Hex_Collector does not need the Mask
						{
							fr = 1; //reuse int fr
							while ( fr <= mskDepth &&
									(MaskSt.getSliceLabel(fr) != null) &&
									!(MaskSt.getSliceLabel(fr).contains(label)) )
									// !( IdataSt.getSliceLabel(fr).startsWith(label) ) && 
									// !( label.startsWith(IdataSt.getSliceLabel(fr))) )
							{	fr++;}
							mask_found = ( (fr <= mskDepth) && (MaskSt.getSliceLabel(fr) != null) );	
							if( (fr <= mskDepth) && (MaskSt.getSliceLabel(fr)==null)) //first test is redundant
							{
								MaskSt.setSliceLabel(label, fr);
								Mask.setStack(MaskSt); //update this one immediately since it may be also be reordered
							}
							
							if( (fr <= mskDepth) && (MaskSt.getSliceLabel(fr) != null) ) //both test are redundant
							{
									
								if ( reorder_Idata )
								{
									//handles to identified frame
									short[] pixels = (short[])MaskSt.getPixels(fr);
									String msklbl = MaskSt.getSliceLabel(fr);
									//put identified frame (fr) to target frNr
									MaskStO.setPixels(pixels, mkNr );
									MaskStO.setSliceLabel(msklbl, mkNr );
								}
								else
								{	mkNr = fr;}
								
								
							}
							else 
							{	
								if(fr >= mskDepth) 
								{	
									IJ.log("WARNING: could not find mask " + label); 
									//reply = false;	
								}
							}
						}	
					}
				}
				else if ( key.equals("slice") )
				{	
					//silently ignore invalid slice entries, but just use them
					//if the label did not match anything
					int tmpfr = ls.nextInt();
					if(label == null)
					{	label = ("sl"+tmpfr);}
					if(!slice_found)
					{	
						frNr = tmpfr;
						slice_found = true;
						if( IdataSt.getSliceLabel(frNr)==null )
						{
							IdataSt.setSliceLabel(label, frNr);
							Idata.setStack(IdataSt);
						}
						
					}
					if( (Mask != null) && (MaskSt != null) && (!mask_found) )
					{	
						mkNr = tmpfr;
						mask_found = true;
						if( MaskSt.getSliceLabel(mkNr)==null )
						{
							MaskSt.setSliceLabel(label, mkNr);
							Mask.setStack(MaskSt);
						}	
					}		
				}
				else if ( key.equals("dirt") )
				{	free_area = 1.0-ls.nextDouble();}
				else if ( key.equals("graphene") )
				{	free_area = ls.nextDouble();}
				else if ( key.equals("peaks") )
				{	
					peaks = ls.nextInt();
					/*
					if(ls.hasNextInt()) peaks = ls.nextInt();
					else
					{
						System.out.println("numpeak cannot be: " + ls.next());
					}*/
				}
				else if ( key.equals("tuning") )
				{	tuning = ls.nextDouble();}
				else if ( key.equals("tilt") )
				{
					double rotation = ls.nextDouble();
					rotation = rotation % (Math.PI/3); //60°
					if(rotation > Math.PI/6)// stay below +30°
					{	rotation -= Math.PI/3;} 
					if(rotation <= -Math.PI/6)
					{	rotation += Math.PI/3;}// stay above -30°
					tilt = rotation; // -30° < tilt <= +30°
				}
				else if ( key.equals("ellA") )
				{	ea   = ls.nextDouble();}
				else if ( key.equals("ellB") )
				{	eb   = ls.nextDouble();}
				else if ( key.equals("ellPhi") )
				{	ephi= ls.nextDouble();}
				else if ( key.equals("hel") )
				{	hexedgelen = ls.nextDouble();} //in principle also redundant
				else if ( key.equals("offX") )
				{	offset_x = ls.nextDouble();}
				else if ( key.equals("offY") )
				{	offset_y = ls.nextDouble();}		
				else if ( key.equals("Nsf") )
				{	subframes = ls.nextInt();}
				else if ( key.equals("offQ") )
				{	offset_hexq = ls.nextInt();}
				else if ( key.equals("offR") )
				{	offset_hexr = ls.nextInt();}
				else if ( key.equals("time") )
				{	elapsed = ls.nextDouble();}
				else if ( key.equals("merit") )
				{	
					merit = ls.nextDouble();
					old_merit = merit;
				}
				else if ( key.equals("hex_merit") )
				{	hex_merit = ls.nextDouble();}
				else if ( key.equals("uc_avg") )
				{	raw_avg = ls.nextDouble();}
				else if ( key.equals("hex_avg") )
				{	hex_avg = ls.nextDouble(); }
				else if ( key.equals("contrast") )
				{	contrast = ls.nextDouble();}
				else if ( key.equals("positionQ") )
				{	positionQ = ls.nextDouble();}
				else if ( key.equals("mirrorQ") )
				{	mirrorQ = ls.nextDouble();}
				else if ( key.equals("hex_low") )
				{	hex_low = ls.nextDouble();}
				else if ( key.equals("hex_high") )
				{	hex_high = ls.nextDouble();}
				else if ( key.equals("hex_shape") )
				{	hex_shape = ls.nextDouble();}
				else if ( key.equals("moment2") )
				{	moment2 = ls.nextDouble();}
				else if ( key.equals("stageX") )
				{	stageX = ls.nextDouble();}
				else if ( key.equals("stageY") )
				{	stageY = ls.nextDouble();}
				else if ( key.equals("stageN") )
				{	stageN = ls.nextInt();}
				else if ( key.equals("state") )
				{	state = ls.nextInt();
					if (state == 2)
					{	state = 0;}
				}
				else if ( key.equals("imgID") )
				{	imgID = ls.nextInt();}
				else if ( key.equals("hex_std") )
				{	ls.next();} //redundant
				else if ( key.equals("uc_std") )
				{	ls.next();} // redundant
				else if ( key.equals("excent") )
				{	ls.next();} //redundant	
				else if ( key.equals("shapeQ") )
				{	ls.next();} //redundant
				else if ( key.equals("sym") )
				{	ls.next();} //redundant
				else if ( key.equals("<none>") )
				{	ls.next();} //well simply ignore this one
				else
				{	
					IJ.log("ERROR: Unknown key/alias while parsing key: \"" + plain_key+"\" alias: \""+ key +"\"" );
					aliases.remove(plain_key); // we know this one sucks so trash it
					get_an_alias(plain_key);
					continue; //retry with the new key	
				}
			}	
			if(eb > ea)
			{
				double tmp = eb;
				eb = ea;
				ea = tmp;
				ephi = (ephi+0.5*Math.PI)%Math.PI;
			}
			ellipseA = ea;
			ellipseB = eb;
			ellipsePhi = Math.atan2(Math.sin(ephi),Math.cos(ephi));
			if(ellipsePhi > Math.PI/2)
			{	ellipsePhi -= Math.PI;}
			if(ellipsePhi <= -Math.PI/2)
			{	ellipsePhi += Math.PI;}
			update_redundant();
			ls.close();
			return true;
		}
		catch(Exception e)
		{
			System.out.println("Error while scanning");
			System.out.println("KEYS\t" + keyline);
			System.out.println("LINE\t" + line);
			System.out.println("FOR\t" + plain_key + "\taka\t" + key); 
			e.printStackTrace();
			return false;
		}	
	}
	
	private void update_redundant()
	{
		hex_std = hex_avg * contrast;
		//hex_std_scaled = hex_std; 
		raw_std = raw_avg * contrast;
		excent  = ellipseA / ellipseB;
		shapeQ  = mirrorQ  * positionQ;
		sym = hex_shape / moment2;
		if( unread_label && ( Double.isNaN(stageX) || Double.isNaN(stageY) || (imgID == -1) ) && 
				(Idata != null) && (label != null) ) //only effective if an entry was not read earlier
		{
			try
			{
				Scanner ls = new Scanner(label);
				Pattern ptrn = Pattern.compile("_");
				ls.useDelimiter(ptrn);
				imgID = ls.nextInt();
				stageX = ls.nextDouble();
				stageY = ls.nextDouble();
				if(ls.hasNextInt())
				{
					stageN = ls.nextInt();
				}
			}
			catch (Exception e)
			{
				if(label_ok)
				{
					System.out.println("Exception while parsing label: \""+ label + "\" of frame: " + frNr);
					System.out.println("suppressing further ill formated labels");
					label_ok = false;
				}
				//e.printStackTrace();
				imgID = frNr;
				stageX = 0;
				stageY = 0;
				stageN = -1;
				
			}
			unread_label = false;
		}
	}
	
	/*
	public HexSampleReport( HexSampleReport other)
	{
		elapsed = other.elapsed;
		hexedgelen = other.hexedgelen;
		rot = other.rot;
		offset_x = other.offset_x;
		offset_y = other.offset_y;
		merit_edge = other.merit_edge;
		merit_xyr = other.merit_xyr;
		modelsize = other.modelsize;	
	}
	*/ 

}
