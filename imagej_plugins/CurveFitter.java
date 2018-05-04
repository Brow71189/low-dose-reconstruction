import java.util.Arrays;
import java.util.Random;


class CurveFitter
{
	Random rand = new Random();
	
	private double[] pts_x = null;
	private double[] pts_y = null;
	//private double[] pts_y2 = null;
	//private double[] pts_swap = null; 
	
	private double const_avg = Double.NaN;
	private double avg_x = Double.NaN;
	private double avg_y = Double.NaN;
	private double lin_d = Double.NaN; 
	private double lin_k = Double.NaN;
	private double o1_1 = Double.NaN;
	private double o1_0 = Double.NaN;
	
	private double root_2 = Double.NaN;
	private double root_1 = Double.NaN;
	private double root_0 = Double.NaN; 
	
	private int pts_len = -1;
	String[] models = {	"const" , "linear", "root" };
	private double[] merits = {Double.NaN, Double.NaN, Double.NaN};
	int mode = 1; //linear
	
	CurveFitter(){}
	
	public void set_pts(final double[] x, final double[] y)
	{
		pts_len = x.length;
		if(pts_len != y.length)
		{
			throw new IllegalStateException("x and y must be of same length");
		}
		pts_x = Arrays.copyOf(x,pts_len);
		pts_y = Arrays.copyOf(y,pts_len);
		update_avgs();	
	}
	
	private void update_avgs()
	{
		double a_bar = 0.0;
		double s_bar = 0.0;
		for(int i=0; i < pts_len; ++i)
		{
			a_bar += pts_x[i];
			s_bar += pts_y[i];
			//pts_y2[i] = Math.pow(pts_y[i],2);
		}
		avg_x = a_bar / pts_len;
		avg_y = s_bar / pts_len;
		//System.out.println("avg x,y: " + avg_x + "," + avg_y);
	}
	
	
	public int get_len()
	{	return pts_len;}
	
	
	public void set_model(final int new_mode)
	{
		if( (new_mode < 0) || (new_mode >= models.length) )
		{
			throw new IllegalStateException("unknown mode: " + new_mode);
		}
		mode = new_mode;	
	}
	
	public void set_model(final String modename)
	{
		int new_mode = 0;
		while( (new_mode < models.length) && (!modename.equals(models[new_mode])) )
		{	++new_mode;}
		set_model(new_mode);
	}
	
	public String get_model()
	{
		return models[mode];
	}
	
	public int get_model_id()
	{
		return mode;
	}
	
	public String[] get_all_models()
	{
		String[] all_models = Arrays.copyOf(models,models.length+1);
		all_models[models.length] = "auto";
		return all_models;
	}
	
	public double get_val(final double x)
	{
		switch(mode)
		{
			case 0:
				return avg_y;
			case 1:
				return o1_1 * (x - avg_x) + o1_0;
			case 2:
				return Math.sqrt( root_2 * Math.pow(x-avg_x,2) + root_1 * (x - avg_x) + root_0 );	
		}
		return Double.NaN;
	}
	
	public double get_val(final int m, final double x)
	{
		switch(m)
		{
			case 0:
				return avg_y;
			case 1:
				return o1_1 * (x - avg_x) + o1_0;
			case 2:
				return Math.sqrt( root_2 * Math.pow(x-avg_x,2) + root_1 * (x - avg_x) + root_0 );	
		}
		return Double.NaN;
	}
	
	private double residual()
	{
		double error = 0.0;
		for(int i = 0; i < pts_len; ++i)
		{	error += Math.pow(pts_y[i]-get_val(pts_x[i]),2); }
		return error;
	}
	
	public double[] get_parameters()
	{
		switch(mode)
		{
			case 0:
				return new double[]{avg_x,avg_y};
			case 1:
				return new double[]{avg_x, o1_0, o1_1};
			case 2:
				return new double[]{avg_x,root_0, root_1, root_2};	
		}
		return null;
	}
	
	public double get_zero()
	{
		return -lin_d/lin_k;
	}
	
	public String get_formula()
	{
		return get_formula(mode);
	}
	
	public String get_formula(int m)
	{
		switch(m)
		{
			case 0:
				return "y=" + avg_y;
			case 1:
				return "y=" + o1_0 + "+" + o1_1 +"*(x-" + avg_x + ")" + "  y=0->x=" + (-lin_d/lin_k);
			case 2:
				return "y=" + "{ " + root_0 + "+" + root_1 + "*(x-" + avg_x + ")" +
				"+" + root_2 + "*(x-" + avg_x + ")^2 }^0.5";	
		}
		return null;
	}
	
	
	
	public double get_merit()
	{
		return merits[mode];
	}
	
	public double[] get_all_merits()
	{
		return merits;
	}

	public double fit(final int new_model)
	{
		set_model(new_model);
		return fit();
	}
	
	public double fit(final String modelname)
	{
		if(modelname.equals("auto"))
		{
			double best_merit = Double.NaN;
			int best_m = -1;
			boolean first = true;
			for(int m = 0; m < models.length; ++m)
			{
				set_model(m);
				fit();
				//System.out.println("merit: " + merits[m]);
				if(first || merits[m] < best_merit)
				{	
					best_merit = merits[m];
					best_m     = m;	
				}
			}
			set_model(best_m);
			return best_merit;
		}
		else
		{
			set_model(modelname);
			return fit();
		}
	}
	
	
	public double fit()
	{
		switch(mode)
		{
			case 0:
			//	fit_const();
			break;
			case 1:
				fit_linear();	
			break;
			case 2:
				fit_root();
			break;
			default:
				throw new IllegalStateException("unknown mode: " + mode);

		}
		return ( merits[mode] = Math.sqrt( residual()/pts_len ) );
	}

	private void fit_linear()
	{
		//double a_bar = avg_x;
		//double s_bar = avg_y;
		
		double xx_var = 0.0; //xx_var
		//double ss_bar = 0.0;
		double xy_var = 0.0; //xy_var
		
		for(int i=0; i < pts_len; ++i)
		{
			xx_var += Math.pow(pts_x[i] - avg_x,2);
			//ss_bar += Math.pow(pts_y[i] - avg_y,2);
			xy_var += (pts_y[i] - avg_y) * (pts_x[i] - avg_x);
		}
		
		lin_k = xy_var/xx_var;
		lin_d = avg_y - lin_k * avg_x;
		o1_1 = lin_k;
		o1_0 = avg_y;
		/*if(lin_k > 0.8)
		{
			System.out.println("WARNING   excessive slope of std vs avg");
			System.out.println("linear y=kx+d\tk: " + lin_k + "\td: " + lin_d);
		}*/
		return;
	}
	
	private void write_root_parameters(double vec[])
	{
		root_0 = vec[0];
		root_1 = vec[1];
		root_2 = vec[2];
		return;
	}
	
	private double[] add(double[] a, double[] b)
	{
		double c[] = new double[a.length];
		for(int i=0; i < a.length; ++i)
		{	c[i] = a[i] + b[i];}
		return c;
	}
	
	private void expand(double[] step )
	{
		for(int i=0; i < step.length; ++i)
		{	step[i] *= 1.35; } //
	}
	
	private void shrink(double[] step )
	{
		for(int i=0; i < step.length; ++i)
		{	step[i] *= -0.65; } //
	}		
	
	private void fit_root()
	{
		fit_linear();
		root_0 = o1_0 * o1_0;
		root_1 = 2 * o1_0 * o1_1;
		root_2 = o1_1 * o1_1;
		double droot_0 = 0.01 * root_0 * (1+rand.nextDouble());
		double droot_1 = 0.01 * root_1 * (1+rand.nextDouble());
		double droot_2 = 0.01 * root_2 * (1+rand.nextDouble());
		if(rand.nextBoolean())
		{
			droot_0 = -droot_0;	
		}
		if(rand.nextBoolean())
		{
			droot_1 = -droot_1;	
		}
		if(rand.nextBoolean())
		{
			droot_2 = -droot_2;	
		}
		
		double[] vec0 = {root_0, root_1, root_2};
		final int pars = vec0.length;
		double[] vec1 = Arrays.copyOf(vec0,pars);
		double res0 = residual();
		final double res_init = res0;
		double res1 = res0;
		double[] step0 = {droot_0, droot_1, droot_2};
		int expands = 0;
		double recent_progress = 0.0;
		int recent_expands = 0;
		/*
		System.out.println("sqrt(a*(x-x0)^2+b*(x-x0)+c)\tres: " +
							res0 + "\n" +
							"\ta: " + vec0[2] + "," + step0[2] +
							"\tb: " + vec0[1] + "," + step0[1] +
							"\tc: " + vec0[0] + "," + step0[0]   );
		*/
		do
		{	
			recent_expands = 0;
			double res_recent = res0;
			for(int trial=0; trial < 60; ++trial)
			{
				int m = trial%pars;
				if(trial > pars * pars)
				{	m = rand.nextInt(pars);}
				double[] step = new double[pars];
				step[m] = step0[m];
				int shrinks = 0;
				while(shrinks + expands < 25)
				{
					vec1 = add(vec0,step);
					write_root_parameters(vec1);
					res1 = residual();
					if( (res1 < res0) && (root_2 > 0) && (root_1 > 0) && (root_0 > 0) )
					{
						vec0 = Arrays.copyOf(vec1,pars);
						res0 = res1;
						++expands;
						expand(step);
					}
					else
					{
						++shrinks;
						shrink(step);
					}
				}	
				step0[m] = step[m];
			}
			recent_progress = res_recent/res0-1.0;
			/*
			System.out.println("sqrt(a*(x-x0)^2+b*(x-x0)+c)\tres,prog: " +
								res0 + "," + recent_progress + "\n" +
							   "\ta: " + vec0[2] + "," + step0[2] +
							   "\tb: " + vec0[1] + "," + step0[1] +
							   "\tc: " + vec0[0] + "," + step0[0]   );
		   */ 
		}
		while( (recent_expands > 2) || (recent_progress > 0.01) );
		write_root_parameters(vec0);
		return;	
	}
	

}
