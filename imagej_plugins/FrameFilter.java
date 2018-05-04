import ij.*;
import java.util.Scanner;
import java.util.Vector;
import java.util.regex.Pattern;
import java.io.File;
import java.io.FileReader;

public class FrameFilter
{
	private class Bounds
	{
		private double lower_bound = Double.NaN;
		private double upper_bound = Double.NaN;
		private boolean invert = false;
		
		public Bounds(double low, double high, boolean inv)
		{
			lower_bound = low;
			upper_bound = high;
			invert = inv;
		};
		
		public boolean is_inside(double val)
		{
			boolean inside = true;
			if(!Double.isNaN(lower_bound))
			{	inside &= lower_bound < val;}
			if(!Double.isNaN(upper_bound))
			{	inside &= upper_bound >= val;}
			return invert^inside;
		}
	}
	
	private Vector< Vector<Bounds> > filters = null;
	private String[] headers = null;
	private int active_filters = 0;
	
	public FrameFilter(File filter_file)
	{
		headers = HexSampleReport.headerline.split("\\s+");
		filters = new Vector< Vector<Bounds> >(headers.length); 
		for(int i=0; i < headers.length; ++i)
		{
			filters.add(i, new Vector<Bounds>() );
		}
		Scanner ffs = null; 
		boolean ffok = (filter_file != null) && (filter_file.isFile() && filter_file.setReadable(true)); 
		if(ffok)
		{	
			try
			{
				ffs = new Scanner(new FileReader(filter_file));
				IJ.log("reading " + filter_file.getPath());
			}
			catch (Exception e)
			{
				ffok = false;
				System.out.println("Error reading: " + filter_file.getPath() );
				e.printStackTrace();
			}	
		}
		if(!ffok)
		{	System.out.println("Could not open: " + filter_file.getPath());}
		else
		{ 
			while(ffs.hasNextLine())
			{
				String key = null;
				double low = Double.NaN;
				double high = Double.NaN;
				boolean inv = false;
				
				String line = ffs.nextLine();
				line = line.trim();
				//skip comments and empty lines
				while ( ffs.hasNextLine() && ( (line.startsWith("#") || line.length() == 0) ) )
				{
					line = ffs.nextLine();
					line = line.trim();				
				}
				Scanner ls = new Scanner(line);
				if(ls.hasNext())
				{
					String plain_key = ls.next().trim();
					if(!plain_key.startsWith("#"))
					{
						inv = plain_key.startsWith("!");
						key = plain_key.substring(inv?1:0);
						//We could get funny redirections here by looping ;)
						if(HexSampleReport.aliases.containsKey(plain_key))
						{	key = HexSampleReport.aliases.get(plain_key);}
						if(ls.hasNextDouble())
						{
							low = ls.nextDouble();
							if(ls.hasNextDouble())
							{	high = ls.nextDouble();}
						}
						int id = 0;
						while( (id < headers.length) && !headers[id].equals(key) ) {++id;}
						if( id < headers.length ) 
						{
							filters.get(id).add( new Bounds(low,high,inv) );
							++active_filters;
						}
						else
						{
							System.out.println("Unrecognized key: " + key);
							System.out.println(line);
						}
					}
				}	
			}
			if(active_filters == 0)
			{
				IJ.log("HINT: in filter files any blank lines or lines starting with # are ignorred");
				IJ.log("HINT: valid keys can be found in frame_continue.txt");
				IJ.log("HINT: any key is followed by an exclusive lower and, if present, an inclusive upper bound");
				IJ.log("HINT: a ! before the key will invert the condition");
				IJ.log("HINT: multiple entries with the same key will be ORed");
				IJ.showMessage("There were not any Filters defined in " + filter_file.getName());
			}
			IJ.log("active filters: " + active_filters);
		}	
	}

	public boolean apply_filter( HexSampleReport report)
	{
		boolean passes = true;
		String[] reps = report.toString().split("\\s+");
		for(int i = 0; i < filters.size(); ++i) 
		{
			if(!filters.get(i).isEmpty()) //this value is filtered
			{
				try
				{
					double val = Double.parseDouble(reps[i]);	
					boolean any_true = false;
					for(int ind = 0; ind < filters.get(i).size() && (!any_true) ; ++ind)
					{	any_true |= filters.get(i).get(ind).is_inside(val);}
					passes &= any_true;
				}
				catch( Exception e )
				{
					System.out.println("failed to parse: \"" + reps[i] + "\" as a value of " +  headers[i]);
					e.printStackTrace();
				}
			}
		}
		return passes;
	}
}
