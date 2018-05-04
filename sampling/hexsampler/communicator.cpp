
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <sys/file.h>
#include <inttypes.h>
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "globals.hpp"
#include "basis.hpp"
#include "ellipse.hpp"
#include "communicator.hpp"
#include "unitcell.hpp"
#include "offset.hpp"
#include "optimizer.hpp"
#include "hex_optimizer.hpp"
#include "hexagons.hpp"



//private headers

int clean_pts( 0 );

double getdouble(); //get a binary Big Endian double from stdin
int getFrame(bool mask_received);
int getMask(void);
int getTarget(void); 
void send_uc_pos(void);
void send_uc_sum(void);
void send_ruc_sum(void);
void send_hexes(void);
void send_sf_sum(void);
void send_sf(void);




double getdouble()  //we have to reverse the byte order from JAVA (at least for linux/gcc)
{
	union 
	{
		double dbl;
		char byts[sizeof(double)/sizeof(char)];  	
	};
	for(size_t i = sizeof(double)/sizeof(char); i > 0u; --i)
	{		byts[i-1] = getchar();}
	return dbl;
}

/*
//This does not work, maybe a byte alignment issue with the binary stream
void senddouble(double val)  //we have to reverse the byte order for JAVA (at least from linux/gcc)
{
  char *ret = (char*)&val; //quick and dirty
  putchar(ret[7]);
  putchar(ret[6]);
  putchar(ret[5]);
  putchar(ret[4]);
  putchar(ret[3]);
  putchar(ret[2]);
  putchar(ret[1]);
  putchar(ret[0]);

}
*/
void senddouble(double val)  //we have to reverse the byte order for JAVA (at least from linux/gcc)
{
	union 
	{
		double dbl;
		char byts[sizeof(double)/sizeof(char)];  	
	};
	dbl = val;
  
	for(size_t i = sizeof(double)/sizeof(char); i > 0u; --i)
	{	putchar(byts[i-1]);}
 
}

void sendshort(unsigned short val)
{
	/*
	static int death_counter(100);
	const unsigned char high( (val & 0xFF00) >> 8 );
	const unsigned char low( (val & 0x00FF) );
	const unsigned short check( 256*high + low );
	const unsigned short chk2((high << 8)|low );
	
	if( (check != val) || (chk2 != val) )
	{
		fprintf(logfile,"%d != 256*%d +%d = %d\n", val, high, low, check);
		fflush(logfile);
		if(--death_counter < 0)
		{	abort(); }	
	}
	*/
	
	putchar( (val & 0xFF00) >> 8 ); //high byte of short value
	putchar( (val & 0x00FF) );      //low byte of short value
	//printf("%d\n",val);  //DEBUG textrepresentation
}



void send_uc_pos(void)
{
	int64_t chk_sum( 0 ); //java long is 64 bit
	long sum( 0 );
	const int *v(uc_pos);
	
	for(int ind( 0 ); ind < uc_area; ++ind)
	{	
		sum += *(v++);
	}
	const double scale( 1000.0/( (double)sum/(double)uc_area) );
	fprintf(logfile, "offering uc_pos\n");fflush(logfile);
	if(!has_console)
	{
		//wait until the master is ready to receive
		fflush(stdout);
		globalCommand->assertCommand("SendPixels", 1);
	}
	
	printf("BeginPixels\t%d\n", uc_area);
	for(int ind(0); ind < uc_area; ++ind)
	{
		const int x( ind % uc_res );
		const int y( ind / uc_res );
		
		const int val( (int) ( scale * ( uc_pos_val(x, y) ) ) );
		sendshort( val );
		//printf("%d\n",val);
		chk_sum += val;
	}
	fflush(stdout);
	printf("EndPixels\t%" PRId64 "\n", chk_sum);
	fprintf(logfile, "uc_pos chk_sum: %" PRId64 "\n", chk_sum);fflush(logfile);
	//assert(chk_sum != 0);
}

void send_uc_sum(void)
{
	int64_t chk_sum( 0 );
	long sum( 0 );
	for(int ind( 0 ); ind < uc_area; ++ind)
	{	sum += uc_sum[ind];}
	const double scale( 1000.0*(double)uc_area/(double)sum );
	fprintf(logfile, "offering uc_sum\n");fflush(logfile);
	if(!has_console)
	{
		//wait until the master is ready to receive
		fflush(stdout);
		globalCommand->assertCommand("SendPixels", 1);
	}
	printf("BeginPixels\t%d\n", uc_area);
	for(int ind( 0 ); ind < uc_area; ++ind)
	{
		const int x( ind % uc_res );
		const int y( ind / uc_res );
		const int val = (int)(scale * (  uc_sum_val(x, y) ));
		sendshort( val );
		//printf("%d\n",val);
		chk_sum += val;
	}
	fflush(stdout);
	printf("EndPixels\t%" PRId64 "\n", chk_sum);
	fprintf(logfile, "uc_sum chk_sum: %" PRId64 "\n", chk_sum);fflush(logfile);
	//assert(chk_sum != 0);
}

void send_uc_mini(void)
{
	int64_t chk_sum( 0 );
	long sum( 0 );
	for(int ind( 0 ); ind < uc_mini_len; ++ind)
	{	sum += uc_mini[ind];}
	const double scale( 0.75*1000.0*(double)uc_mini_len/(double)sum );
	fprintf(logfile, "offering uc_mini\n");fflush(logfile);
	if(!has_console)
	{
		//wait until the master is ready to receive
		fflush(stdout);
		globalCommand->assertCommand("SendPixels", 1);
	}
	printf("BeginPixels\t%d\n", uc_mini_len);
	for(int ind( 0 ); ind < uc_mini_len; ++ind)
	{
		const int val = (int)(scale * uc_mini[ind] );
		sendshort( val );
		//printf("%d\n",val);
		chk_sum += val;
	}
	fflush(stdout);
	printf("EndPixels\t%" PRId64 "\n", chk_sum);
	fprintf(logfile, "uc_mini chk_sum: %" PRId64 "\n", chk_sum);fflush(logfile);
}

void send_ruc_sum(void)
{
	int64_t chk_sum( 0 );
	long sum( 0 );
	for(int ind( 0 ); ind < uc_area; ++ind)
	{	sum += ruc_sum[ind];}
	const double scale( 1000.0*(double)uc_area/(double)sum );
	fprintf(logfile, "offering ruc_sum\n");fflush(logfile);
	if(!has_console)
	{
		//wait until the master is ready to receive
		fflush(stdout);
		globalCommand->assertCommand("SendPixels", 1);
	}
	
	printf("BeginPixels\t%d\n", uc_area);
	for(int ind( 0 ); ind < uc_area; ++ind)
	{
		const int x( ind % uc_res );
		const int y( ind / uc_res );
		const int val = (int)( scale * ( ruc_sum_val(x, y) ) );
		
		sendshort( val );
		//printf("%d\n",val);
		chk_sum += val;
	}
	fflush(stdout);
	printf("EndPixels\t%" PRId64 "\n", chk_sum);
	fprintf(logfile, "ruc_sum chk_sum: %" PRId64 "\n", chk_sum);fflush(logfile);
	//assert(chk_sum != 0);
}



void send_hexes(void)
{
	int64_t chk_sum( 0 );
	fprintf(logfile, "offering hexframe\n");fflush(logfile);
	if(!has_console)
	{
		//wait until the master is ready to receive
		fflush(stdout);
		globalCommand->assertCommand("SendPixels", 1);
	}
	
	printf("BeginPixels\t%d\n", hexframelen);
	const int *hex( hexarea );
	//const int *hex(hex_mask);//DEBUG send hex_mask instead of hexarea 
	for(int ind = 0; ind < hexframelen; ++ind)
	{
		int val( *(hex++) );
		if( val == -1 )
		{	val = 0;}
		else if( val == dirt_val)
		{	val = 1;} //usually invisible unless contrast is set to 0..1
		sendshort( val );
		//printf("%d\n",val);
		chk_sum += val;
		//Maybe that is the reason why the clusterhead crashes???
		/*
		if( (1+ind)%4096 == 0 )
		{
			fflush(stdout);
			std::this_thread::sleep_for( std::chrono::milliseconds(2) );
		}
		*/ 
	}
	fflush(stdout);
	printf("EndPixels\t%" PRId64 "\n", chk_sum);
	fprintf(logfile, "hexframe chk_sum: %" PRId64 "\n", chk_sum);fflush(logfile);
	//assert(chk_sum != 0);
}

void send_sf_sum(void)
{
	const int modelarea( modelsize*modelsize );
	const bool do_scaling ( (op_mode == 0) || (op_mode == 3) );
	int64_t chk_sum( 0 );
	long sum( 0 );
	for(int ind( 0 ); ind < modelarea; ++ind)
	{	sum += sf_sum[ind];}
	const double scale( do_scaling ? ( (1000.0/0.75)*(double)modelarea/(double)sum ) : 1.0 );
	fprintf(logfile, "offering sf_sum\n");fflush(logfile);
	if(!has_console)
	{
		//wait until the master is ready to receive
		fflush(stdout);
		globalCommand->assertCommand("SendPixels", 1);
	}
	printf("BeginPixels\t%d\n", modelarea);
	const int *sfs(sf_sum);
	for(int ind(0); ind < modelarea; ++ind)
	{
		//const int val( ( *(sfs++)+(Num_sf/2) )/Num_sf ); //rounded int average
		const int val( do_scaling ? (int)(*(sfs++)*scale) : ( ( *(sfs++)+(Num_sf/2) )/Num_sf ) );
		sendshort( val );
		//printf("%d\n",val);
		chk_sum += val;
	}
	fflush(stdout);
	printf("EndPixels\t%" PRId64 "\n", chk_sum);
	fprintf(logfile, "sf_sum chk_sum: %" PRId64 "\n", chk_sum);fflush(logfile);
	//assert(chk_sum != 0);
}

void send_sf(void)
{
	const int modelarea( modelsize*modelsize );
	fprintf(logfile, "offering subframes:%d\n", Num_sf);fflush(logfile);	
	for(int sfr(0); sfr < Num_sf; ++sfr)
	{
		int64_t chk_sum( 0 );
		if(!has_console)
		{
			//wait until the master is ready to receive
			fflush(stdout);
			globalCommand->assertCommand("SendPixels", 1);
		}
		printf("BeginPixels\t%d\n", modelarea);
		const int *sfs(subframes[sfr]);
		for(int ind(0); ind < modelarea; ++ind)
		{
			const int val( *(sfs++) );
			sendshort( val );
			//printf("%d\n",val);
			chk_sum += val;
		}
		fflush(stdout);
		printf("EndPixels\t%" PRId64 "\n", chk_sum);
		//fprintf(logfile, ".");fflush(logfile);
		//assert(chk_sum != 0);
	}
	fprintf(logfile, "sent sfubframes: %d\n", Num_sf);fflush(logfile);	
}

int getFrame(bool mask_received)
{
	framelen = impWidth * impHeight;
	assert(framelen > 0);
	delete[] frame;
	frame = new int[framelen];
	memset(frame,0,framelen*sizeof(int));
	if(!mask_received)
	{
		delete[] mask;
		mask = new unsigned char[framelen];
		memset(mask,0,framelen); //all flag_hole by default
	}
	fprintf(logfile, "awaiting Idata\n");fflush(logfile);
	const char *head;
	const char *args;
	globalCommand->next();
	head = globalCommand->getHead().c_str();
	args = globalCommand->getArgs().c_str();
	int checksum( 0 );
	long sum( 0 );
	int arg( -1 );
	if(!mask_received)
	{
		clean_X = 0.0;
		clean_Y = 0.0;
	}
	double iW_2 (impWidth/2);
	double iH_2 (impHeight/2);
	bool arg_ok( sscanf( args, "%d", &arg) );
	int *fr(frame);
	if (arg_ok && (strcmp(head, "BeginBinary") == 0) && (arg == 2*framelen) )
	{	
		//globalCommand->assertCommand("BeginBinary", 2 * framelen );
		
		for (int ind = 0; ind < framelen; ++ind)
		{
			unsigned short val = getchar() << 8;
			val |= getchar();
			checksum += val;
			if(!mask_received || (mask[ind] == 1) )
			{	sum += val;}
			*(fr++) = val;		
			if( (!mask_received) && (val < dirt_val) )
			{
				++clean_pts;
				const double x( ind % impWidth - iW_2);
				const double y( ind / impWidth - iH_2);
				clean_X += x;
				clean_Y += y;
				mask[ind] = 1; // 1 .. flag_graphene
			}
		}
		globalCommand->assertCommand("EndBinary", checksum );
		//globalCommand->next();
    }
    else if(arg_ok && (strcmp(head, "BeginPixels") == 0) && (arg == framelen) )
    {
		//int *fr(frame);
		bool first = true;
		for (int ind = 0; ind < framelen; ++ind)
		{
			unsigned short val( 0 );
			if( (scanf("%hu\n",&val) != 1) && first)
			{
				printf("ERROR in getFrame() while reading text pixel nr %d",ind);
				fflush(stdout);
				first = false;
			}
			checksum += val;
			sum += val;
			if(!mask_received || (mask[ind] == 1) )
			{	sum += val;}
			*(fr++) = val;	
			if( ( (!mask_received) && (val < dirt_val) ) )
			{
				++clean_pts;
				const double x( ind % impWidth - iW_2);
				const double y( ind / impWidth - iH_2);
				clean_X += x;
				clean_Y += y;
				mask[ind] = 1;
			}	
		}
		globalCommand->assertCommand("EndPixels", checksum );
	}
    else
    {	
		printf("ERROR in getFrame() expected either BeginBinary(%d) or BeginPixels(%d)\n", 2 * framelen, framelen);
		fflush(stdout);
		fprintf(logfile,"ERROR in getFrame() expected either BeginBinary(%d) or BeginPixels(%d)\n", 2 * framelen, framelen);
		fflush(stdout);
		exit(0);	
	}
    if(!mask_received)
    {
		if(clean_pts > 0)
		{
			clean_X = clean_X/(double)clean_pts + (double)iW_2;
			clean_Y = clean_Y/(double)clean_pts + (double)iH_2;
		}
		else
		{
			clean_X = iW_2;
			clean_Y = iH_2;
		}
	}
    //check_memory();
    //assert(MCHECK_OK == mprobe(mask));
    //assert(MCHECK_OK == mprobe(frame));
    fprintf(logfile, "Idata chk_sum: %d\n", checksum);fflush(logfile);
    fprintf(logfile, "Idata sum: %ld, avg: %lf\n", sum, ((double)sum)/(mask_received?clean_pts:framelen) );
    return checksum;
}

int getMask(void)
{
	const int mskWidth(impWidth/mask_scaling);
	const int mskHeight(impHeight/mask_scaling);
	
	const int masklen( mskWidth * mskHeight);
	assert(masklen > 0);
	delete[] mask;
	mask = new unsigned char[masklen];
	memset(mask,0,masklen);  //all flag_hole by default
	fprintf(logfile, "awaiting Mask\n");fflush(logfile);
	const char *head;
	const char *args;
	globalCommand->next();
	head = globalCommand->getHead().c_str();
	args = globalCommand->getArgs().c_str();
	int checksum( 0 );
	long sum( 0 );
	int arg( -1 );
	clean_pts = 0;
	clean_X = 0.0;
	clean_Y = 0.0;
	double iW_2 (mskWidth/2);
	double iH_2 (mskHeight/2);
	bool arg_ok( sscanf( args, "%d", &arg) );
	unsigned char *fr(mask);
	if (arg_ok && (strcmp(head, "BeginBinary") == 0) && (arg == masklen) )
	{	
		//globalCommand->assertCommand("BeginBinary", 2 * framelen );
		
		for (int ind = 0; ind < masklen; ++ind)
		{
			unsigned char val = (unsigned char)getchar();
			checksum += val;
			sum += val;
			*(fr++) = val;		
			if(val == 1)
			{
				++clean_pts;
				const double x( ind % mskWidth - iW_2);
				const double y( ind / mskWidth - iH_2);
				clean_X += x;
				clean_Y += y;
			}
		}
		globalCommand->assertCommand("EndBinary", checksum );
		//globalCommand->next();
    }
    else if(arg_ok && (strcmp(head, "BeginPixels") == 0) && (arg == masklen) )
    {
		//int *fr(frame);
		bool first = true;
		for (int ind = 0; ind < masklen; ++ind)
		{
			unsigned char val( 0 );
			if( (scanf("%hhu\n",&val) != 1) && first)
			{
				printf("ERROR in getMask() while reading text pixel nr %d",ind);
				fflush(stdout);
				first = false;
			}
			checksum += val;
			sum += val;
			*(fr++) = val;
			if(val == 1)
			{
				++clean_pts;
				const double x( ind % mskWidth - iW_2);
				const double y( ind / mskWidth - iH_2);
				clean_X += x;
				clean_Y += y;
			}	
		}
		globalCommand->assertCommand("EndPixels", checksum );
	}
    else
    {	
		printf("ERROR in getFrame() expected either BeginBinary(%d) or BeginPixels(%d)\n", masklen, masklen);
		fflush(stdout);
		fprintf(logfile,"ERROR in getFrame() expected either BeginBinary(%d) or BeginPixels(%d)\n", masklen, masklen);
		fflush(logfile);
		exit(0);	
	}
    if(clean_pts > 0)
    {
		clean_X = mask_scaling * (clean_X/(double)clean_pts + (double)iW_2);
		clean_Y = mask_scaling * (clean_Y/(double)clean_pts + (double)iH_2);
    }
    else
    {
		clean_X = mask_scaling * iW_2;
		clean_Y = mask_scaling * iH_2;
	}
    //check_memory();
    //assert(MCHECK_OK == mprobe(mask));
    fprintf(logfile, "Mask chk_sum: %d\n", checksum);fflush(logfile);
    return checksum;
}






int getTarget(const bool fake_read)
{
	const int old_uc_mini_len(uc_mini_len); //should always be -1 here
	int * old_bondlength = bondlength;
	//assert(uc_mini_len == -1); //This one passes
	const int len (4*bondlength_Mini*bondlength_Mini); //uc_mini_len is not yet defined berfore uc_mini exists
	uc_mini_len = len; //We need to fake that for using stats from hexagons.hpp
	bondlength = &bondlength_Mini;
	int *uctf( new int[len] );
	memset(uctf,0,len*sizeof(int));
	fprintf(logfile, "awaiting target\n");fflush(logfile);
	globalCommand->next();
	const char *head( globalCommand->getHead().c_str() );
	const char *args( globalCommand->getArgs().c_str() );
	int checksum( 0 );
	int arg( -1 );
	const bool arg_ok( sscanf( args, "%d", &arg)==1 );
	if( arg_ok && (strcmp(head, "BeginBinary") == 0) && (arg == 2*len) )
	{	
		int *uct(uctf);
		for (int ind (0); ind < len; ++ind)
		{
			unsigned short val = getchar() << 8;
			val |= getchar();
			checksum += val;
			*(uct++) = val;
		}
		globalCommand->assertCommand("EndBinary", checksum );
    }
    else if( arg_ok && (strcmp(head, "BeginPixels") == 0) && (arg == len) )
    {
		int *uct(uctf);
		bool first = true;
		for (int ind = 0; ind < len; ++ind)
		{
			unsigned short val( 0 );
			if( (scanf("%hu\n",&val) != 1) && first)
			{
				printf("ERROR in getTarget() while reading pixel nr %d",ind);
				fflush(stdout);
				first = false;
			}
			checksum += val;
			*(uct++) = val;
		}
		globalCommand->assertCommand("EndPixels", checksum );
	}
    else
    {	
		printf("ERROR in getTarget() expected either BeginBinary[%d] or BeginPixels[%d]\n", 2 * len, len);
		fflush(stdout);
		fprintf(logfile, "ERROR in getTarget() expected either BeginBinary[%d] or BeginPixels[%d]\n", 2 * len, len);
		fflush(logfile);
		exit(0);		
	}
    //fprintf(logfile, "* ");fflush(logfile);
    delete[] uc_target;
    uc_target = new double[len];
    //const double iavg ( (checksum > 0) ? (0.75*len)/(double)checksum : 1.0 );
    //fprintf(logfile, "target iavg: %lf ", iavg);fflush(logfile);
    memset(uc_target,0,len*sizeof(double));
    for(int i(0); i < len; ++i)
    {
		uc_target[i] = ((double)uctf[i]); //*iavg //zeros at borders dont matter
	}
	double n_std;
	double n_avg( get_hex_stats( uc_target, n_std) );
	//fprintf(logfile, "new_avg: %lf new_std: %lf \n",  n_avg, n_std);fflush(logfile);
	//normalize_hex( uc_target, n_avg, n_std);
	//n_avg = get_hex_stats( uc_target, n_std);
	fprintf(logfile, "norm_avg: %lf norm_std: %lf \n", n_avg, n_std);fflush(logfile);
	//fprintf(logfile, "* ");fflush(logfile);
    if(fake_read) //we had to do a dake read, but wont actually use the target
    {
		fprintf(logfile, "dismissing uc_target\n");
		delete[] uc_target;
		uc_target = nullptr;
	
	}
    //assert(MCHECK_OK == mprobe(uctf));
    delete[] uctf;
    //assert(MCHECK_OK == mprobe(uc_target));
    //check_memory();
    fprintf(logfile, "checksum: %d\n", checksum);fflush(logfile);
    uc_mini_len = old_uc_mini_len;
    bondlength = old_bondlength;
    return checksum;
}



int read_input(void) //0 .. launch, -1 .. quit now
{
	fprintf(logfile,"receiving input\n"); fflush(logfile);
	int result = 0;
	solidsize = -1;
	const char *head;
	const char *args;
	double new_hel = -1.0;
	double new_tilt = tilt;
	double new_ea = ea;
	double new_eb = eb;
	double new_phi = phi;
	double hexImgScale = 1.0;
	bool mini_len_received = false;
	bool impH_received = false;
	bool impW_received = false;
	bool mask_received = false;
	double dx = offset_X;
	double dy = offset_Y;
	do
	{
		fprintf(logfile,"->\t"); fflush(logfile);
		globalCommand->next();
		head = globalCommand->getHead().c_str();
		args = globalCommand->getArgs().c_str();
		//fprintf(logfile,"at read input \t%s[%s]\n", head, args); fflush(logfile);
		//need to check for "Run" BEFORE getFrame(), getMask() or getTarget() invalidate head and args
		if ( strcmp( head, "Run") == 0 )
		{	
			int val;
			sscanf( args, "%d", &val);
			if (val == 0)
			{	result = -1;}
			break;
		} 
		
		if (strcmp(head, "ThreadNum") == 0)
		{
			sscanf( args, "%d", &threadnum);
		}
		else if (strcmp(head, "Frame") == 0 )
		{
			sscanf( args, "%d", &frameNr);
		}
		else if (strcmp(head, "ImpWidth") == 0 )
		{
			sscanf( args, "%d", &impWidth);
			impW_received = true;
		}
		else if (strcmp(head, "ImpHeight") == 0 )
		{
			sscanf( args, "%d", &impHeight);
			impH_received = true;
		}
		else if (strcmp(head, "MaskScale") == 0 )
		{
			sscanf( args, "%d", &mask_scaling);
		}
		else if (strcmp(head, "HexImgScale") == 0 )
		{
			sscanf( args, "%lf", &hexImgScale);	
		}
		else if (strcmp(head, "ModelSize") == 0 )
		{
			sscanf( args, "%d", &modelsize);
		}
		else if (strcmp(head, "SolidSize") == 0 )
		{
			sscanf( args, "%d", &solidsize);
		}
		else if (strcmp(head, "Resolution") == 0 )
		{
			sscanf( args, "%d", &uc_res);
		}
		else if (strcmp(head, "BondLength") == 0 )
		{
			sscanf( args, "%d", &bondlength_CC);
		}
		else if (strcmp(head, "MiniLength") == 0 )
		{
			mini_len_received  = true;
			sscanf( args, "%d", &bondlength_Mini);
		}
		else if (strcmp(head, "FieldofView") == 0 )
		{
			sscanf( args, "%lf", &fov);
		}
		else if (strcmp(head, "HexEdgeLen") == 0 )
		{
			sscanf( args, "%lf", &new_hel);	
		}
		else if (strcmp(head, "Tilt") == 0 )
		{
			sscanf( args, "%lf", &new_tilt);
		}
		else if (strcmp(head, "EllipseA") == 0 )
		{
			sscanf( args, "%lf", &new_ea);
		}
		else if (strcmp(head, "EllipseB") == 0 )
		{
			sscanf( args, "%lf", &new_eb);
		}
		else if (strcmp(head, "EllipsePhi") == 0 )
		{
			sscanf( args, "%lf", &new_phi);
		}
		else if ( (strcmp(head, "Offset_x") == 0) || (strcmp(head, "OffsetX") == 0) )
		{
			sscanf( args, "%lf", &dx);
		}
		else if ( (strcmp(head, "Offset_y") == 0) || (strcmp(head, "OffsetY") == 0) )
		{
			sscanf( args, "%lf", &dy);
		}
		else if (strcmp(head, "OffsetHexQ") == 0 )
		{
			sscanf( args, "%d", &offset_Q);
		}
		else if (strcmp(head, "OffsetHexR") == 0 )
		{
			sscanf( args, "%d", &offset_R);
		}
		else if (strcmp(head, "NumSF") == 0 )
		{
			sscanf( args, "%d", &Num_sf);
		}
		else if (strcmp(head, "Merit") == 0 )
		{
			sscanf( args, "%lf", &merit);
		}
		else if (strcmp(head, "Hex_Merit") == 0 )
		{
			sscanf( args, "%lf", &hexagonal_merit);
		}
		else if (strcmp(head, "Hex_Shape") == 0 )
		{
			sscanf( args, "%lf", &hex_shape);
		}
		else if (strcmp(head, "Raw_Avg") == 0 )
		{
			sscanf( args, "%lf", &data_avg);
		}
		else if (strcmp(head, "Hex_Avg") == 0 )
		{
			sscanf( args, "%lf", &hex_avg);
		}
		else if (strcmp(head, "Hex_Low") == 0 )
		{
			sscanf( args, "%lf", &hex_low);
		}
		else if (strcmp(head, "Hex_High") == 0 )
		{
			sscanf( args, "%lf", &hex_high);
		}
		else if (strcmp(head, "Contrast") == 0 )
		{
			sscanf( args, "%lf", &contrast);
		}
		else if (strcmp(head, "MinHexContrast") == 0 )
		{
			sscanf( args, "%lf", &min_hex_contrast);
		}
		else if (strcmp(head, "Mirror_Q") == 0 )
		{
			sscanf( args, "%lf", &mirror_quality);
		}
		else if (strcmp(head, "Position_Q") == 0 )
		{
			sscanf( args, "%lf", &position_quality);
		}
		else if (strcmp(head, "Min_Light_Dirt") == 0 )
		{
			sscanf( args, "%lf", &min_light_dirt);
		}
		else if (strcmp(head, "Max_Light_Dirt") == 0 )
		{
			sscanf( args, "%lf", &max_light_dirt);
		}
		else if (strcmp(head, "Light_Dirt_Value") == 0 )
		{
			sscanf( args, "%d", &light_dirt_val);
		}
		else if (strcmp(head, "Retune_Origin") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			retune_origin = (val == 1);
		}
		else if (strcmp(head, "Locking") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			locking = (val == 1);
		}
		else if (strcmp(head, "Skip_Grid") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			skip_grid = (val == 1);
		}
		else if (strcmp(head, "Create_Defect") == 0 )
		{
			sscanf( args, "%d", &defect_type);	
		}
		else if (strcmp(head, "OffsetHexR") == 0 )
		{
			sscanf( args, "%d", &offset_R);
		}
		else if (strcmp(head, "Smoothing") == 0 )
		{
			sscanf( args, "%d", &smooth_passes);
		}
		else if (strcmp(head, "Stability") == 0 )
		{
			sscanf( args, "%d", &init_stability);
		}
		else if (strcmp(head, "PeakPoints") == 0 )
		{
			sscanf( args, "%zu", &top_points);
		}
		else if (strcmp(head, "ReportUC") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			report_uc = (val == 1);
		}
		else if (strcmp(head, "ReportHexImg") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			report_hexes = (val == 1);
		}
		else if (strcmp(head, "ReportSF") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			report_sf = (val == 1);
		}
		else if (strcmp(head, "MarkSF") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			highlight_sf = (val == 1);
		}
		else if (strcmp(head, "Reporting") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			reporting = (val == 1);
		}
		else if (strcmp(head, "HasConsole") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			has_console = (val == 1);
		}
		else if (strcmp(head, "UseMirror") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			use_mirror = (val == 1);
		}
		else if (strcmp(head, "UsePosition") == 0 )
		{
			int val;
			sscanf( args, "%d", &val);
			use_position = (val == 1);
		}
		/*else if (strcmp(head, "RulingMerit") == 0 )
		{
			//ignore this one we have auto
			//sscanf( args, "%d", &ruling_merit);
		}*/
		else if (strcmp(head, "RotateTilt") == 0 )
		{
			sscanf( args, "%d", &rotate_tilt);
		}
		else if (strcmp(head, "Task") == 0 )
		{
			sscanf( args, "%d", &op_mode);
		}
		else if ( strcmp(head, "updateFrameData") == 0)
		{
			if(!(impW_received && impH_received))
			{
				fprintf(stderr,"ERROR a frame cannot be read with unspecified ImpWidth and ImpHeight\n");
				fflush(stdout);
				fprintf(logfile,"ERROR a frame cannot be read with unspecified ImpWidth and ImpHeight\n");
				fclose(logfile);
				exit(0);
			}
			getFrame(mask_received);  //head and args are invalidated
			
		}
		else if ( strcmp(head, "updateMask") == 0)
		{
			if(!(impW_received && impH_received))
			{
				fprintf(stderr,"ERROR a mask cannot be read with unspecified ImpWidth and ImpHeight\n");
				fflush(stdout);
				fprintf(logfile,"ERROR a mask cannot be read with unspecified ImpWidth and ImpHeight\n");
				fclose(logfile);
				exit(0);
			}
			getMask();  //head and args are invalidated
			mask_received = true;
			
		}
		else if ( strcmp(head, "updateUCTarget") == 0)
		{
			if(!mini_len_received)
			{
				fprintf(stderr,"ERROR a target cannot be read with unspecified MiniLength\n");
				fflush(stdout);
				fprintf(logfile,"ERROR a target cannot be read with unspecified MiniLength\n");
				fclose(logfile);
				exit(0);
			}
			getTarget( atoi(args)!=1 );  //head and args are invalidated
			
		}
		else
		{	
			
			fprintf(logfile,"WARNING ignorring unknown initialization command %s(%s)\n", head, args );
			fflush(logfile);
			printf("Message WARNING ignorring unknown initialization command %s(%s)\n", head, args );
			fflush(stdout);	
		}
	} while ( true );
	
	if(result == 0) //no need for that if we would quit anyways
	{
		fprintf(logfile,"processing input ... "); fflush(logfile);
		if (locking && (fd_lock == -1) )
		{
			file_lock = fopen("lock.txt","a");
			if(file_lock != nullptr)
			{
				fd_lock = fileno_unlocked(file_lock);	
				fprintf(logfile,"\nlocking established\n"); fflush(logfile);
				
			}
			else
			{
				locking = false;
				fprintf(logfile,"\nlocking failed\n"); fflush(logfile);
			}
		}
		
		if(!locking)
		{
			fd_lock = -1;
			if(file_lock != nullptr)
			{
				fclose(file_lock);
				file_lock = nullptr; 
			}	
		}
		
		if(solidsize < 1)
		{	solidsize = modelsize;}
		
		if(!mini_len_received)
		{	bondlength_Mini = 3*bondlength_CC;}
		
		if(new_eb > new_ea)
		{
			double old_ea = new_ea;
			new_ea = new_eb;
			new_eb = old_ea;
			new_phi = fmod( phi + M_PI_2l, M_PIl);	
		}
		
		if(new_hel <= 0.0)
		{	new_hel = 2.0 * impWidth/(3.0 * new_ea * (bondlength_CC) * sqrt(3.0) );}
		hel = new_hel;
		tilt = new_tilt;
		offset_X = dx;
		offset_Y = dy;
		offset_XS = offset_X;
		offset_YS = offset_Y;
		
		set_ellipse(new_ea/new_eb, new_phi); //relies on current hel and sets distorted basis
		//check_memory();
		hexImpWidth = (int)(hexImgScale * impWidth);
		hexImpHeight = (int)(hexImgScale * impHeight);
		
		//ensure we have even image dimensions, simplifies odd q,r > cube conversion
		if(hexImpWidth % 2 == 1) { ++hexImpWidth; }
		if(hexImpHeight % 2 == 1) { ++hexImpHeight; }
		//max_coh_len = 7*(hexImpWidth+hexImpHeight)/20; //that will suffice for any reasonable CleanX,CleanY
		//coh_len = max_coh_len;
		init_small_hexes();
		fprintf(logfile,"done\n"); fflush(logfile);
		clean_pts = 0;//reset counter for the next frame
	}
	//check_memory();
	
	return result;
}

void report(void)
{
	if(config.fresh_config)
	{
		bool got_lock(false);
		clock_t now = clock();
		elapsed = double(now-timer)/CLOCKS_PER_SEC;
		if(locking && (fd_lock != -1) )
		{
			flock(fd_lock,LOCK_EX);
			got_lock = true;
		}
		report_results();
		if(got_lock)
		{
			flock(fd_lock,LOCK_UN);
		}	
	}
}

void report_results(void)
{
	++report_num;
	//bool get_new_target( false );
	fprintf(logfile,"offering REPORT on frame %d\n", frameNr); fflush(logfile);
	const char *head;
	const char *args;
	
	if(!has_console)
	{
		bool go( false );
		printf("REPORT\n");
		fflush(stdout);
		do
		{
			globalCommand->next();
			head = globalCommand->getHead().c_str();
			args = globalCommand->getArgs().c_str();
			if(strcmp(head,"Report")==0)
			{
				int val = 0;
				if( (sscanf(args,"%d",&val) == 1) && (val == 1) )
				{
					go = true;
					//fprintf(logfile, " commencing\n"); fflush(logfile);
				}
				else
				{
					fprintf(logfile, " declined\n"); fflush(logfile);
					return; //report is declined
				}
			}
			/*else if(strcmp(head,"updateUCTarget")==0)
			{	
				int val = 0;
				if( (sscanf(args,"%d",&val) == 1) && (val == 1) )
				{
					getTarget();
					if(!has_console)
					{	//actually default behaviour with console
						if(!reporting)
						{	printf("Message "); }
						printf("WARNING encountered delayed target unitcell in stdin\n");
						fflush(stdout);
					}
				}	
			}*/
			else if (strcmp(head,"Run"))
			{
				int val = 0;
				if( (sscanf(args,"%d",&val) == 1) && (val == 0) )
				{
					if(!reporting)
					{	printf("Message "); }
					printf("received Run(0) instead of Report(1)\n");
					fprintf(logfile, "\nreceived Run(0) instead of Report(1)\n");
					busy = false;
					fflush(stdout);
					fclose(logfile);
					//std::remove(logname.c_str());
					exit(0);
				}
			}
		}while(!go);	
	}
	
	printf("Frame\t\t%d\n",		( ( busy && (!has_console) ) ? -1 : frameNr) ); //master checks now if this matches the submitted frame
	//Always send those since their might have been a not yet reported updated before the finalization request
	if(reporting && has_console && (ruling_merit == 0) )
	{	printf("Merit\t\t%lf\t%lf\n",	merit, ( !( sampled_merits.empty() || top5.empty() ) )? sampled_merits[ top5.back() ] : 0.0);}
	else
	{	printf("Merit\t\t%lf\n",	merit);}
	if(reporting && has_console && (ruling_merit == 1) )
	{	printf("Hex_Merit\t%lf\t%lf\n",	hexagonal_merit, ( !( sampled_merits.empty() || top5.empty() ) ) ? sampled_merits[ top5.back() ] : 0.0);}
	else
	{	printf("Hex_Merit\t%lf\n",  hexagonal_merit);}
	
	printf("Contrast\t%lf\n",   contrast);
	printf("Position_Q\t%lf\n", position_quality);
	printf("Mirror_Q\t%lf\n",   mirror_quality);
	printf("Hex_Shape\t%lf\n",    hex_shape);
	printf("Moment2\t\t%lf\n", moment2);
	printf("HexEdgeLen\t%lf\n",	hel);
	printf("Tilt\t\t%lf\n",		tilt);	
	printf("EllipseA\t%lf\n",	ea);
	printf("EllipseB\t%lf\n",	eb);
	printf("EllipsePhi\t%lf\n",	phi);
	printf("UC_Avg\t\t%lf\n",  uc_avg);
	printf("Hex_Low\t\t%lf\n",  hex_low);
	printf("Hex_Avg\t\t%lf\n",  hex_avg);
	printf("Hex_High\t%lf\n",   hex_high);
	
	if(reporting && has_console)
	{
		printf("OffsetX\t\t%lf\t%lf\n",	offset_X, 	offset_XS );
		printf("OffsetY\t\t%lf\t%lf\n",	offset_Y,   offset_YS );
		//printf("next a1,a2\t %lf\t%lf\n", next_offset_a1, next_offset_a2); //no longer used
		if(!p1_and_p2_applied) //they should never ever show up
		{	printf("next p1,p2\t %lf\t%lf\n", next_offset_p1, next_offset_p2);}
	}
	else
	{
		printf("OffsetX\t\t%lf\n",	offset_X);
		printf("OffsetY\t\t%lf\n",	offset_Y);
	}
	
	if(!busy) //these do not change during merit optimization
	{	
		printf("SubFrames\t%d\n",	Num_sf);
		printf("OffsetHexQ\t%d\n",	offset_Q);
		printf("OffsetHexR\t%d\n",	offset_R);
	}
	
	if(report_uc && (uc_len > 0) ) //there ought to be at least one sampled unitcell
	{
		printf("sumUC\t\t%d\n",       uc_area);
		send_uc_sum();
		printf("symUC\t\t%d\n",       uc_area);
		if(busy)
		{	send_uc_pos(); }
		else
		{	send_ruc_sum();}	 
		if( (uc_mini_len > 0) ) 
		{
			if (uc_mini != nullptr)
			{
				printf("miniUC\t\t%d\n",	uc_mini_len);
				send_uc_mini();
			}
			/*if( (!has_console) && busy && (report_num > 1) &&
			( (op_mode == 0) || (op_mode == 3) ) )
			{
				printf("requestingTarget\n");
				get_new_target = true;
			}*/
		}	
	}
	if(!busy && (hexframelen > 0) ) //dont send these during optimization or if they dont exist at all
	{
		const int modelarea( modelsize*modelsize );
		if(report_hexes && (Num_sf > 0) && (sf_sum != nullptr) )
		{
			printf("sfSum\t%d\n",	modelarea);
			send_sf_sum();
		}
		if(report_hexes && (Num_sf > 0) && (hexarea != nullptr)  )
		{
			printf("hexImg\t%d\n",	hexframelen);
			send_hexes();
		}
		if(report_sf && (Num_sf > 0) && (subframes != nullptr) )
		{
			printf("sfStack\t%d\n", modelarea);
			send_sf();
		}
	}
	printf("Elapsed\t\t%lf\n",		elapsed); //elapsed also signals the end of the report, it must be sent last
	fflush(stdout);
	config.fresh_config = false;
	fprintf(logfile,"Report submitted\n"); fflush(logfile);
	
	/*if(get_new_target)
	{
		fprintf(logfile,"receiving new target ... "); fflush(logfile);
		bool done (false);
		do
		{
			globalCommand->next();
			head = globalCommand->getHead().c_str();
			args = globalCommand->getArgs().c_str();
			if(strcmp(head,"updateUCTarget")==0)
			{	
				done = true;
				int val = 0;
				if( (sscanf(args,"%d",&val) == 1) && (val == 1) )
				{
					getTarget();
					if(reporting)
					{	printf("updated target unitcell\n");fflush(stdout);}
				}	
			}
			else if(reporting)
			{	
				printf("Error expected updateUCTarget(0/1) instead of: %s[%s]\n",head,args);
				fflush(stdout);
				usleep(50000);
			}
		}
		while (!done);	
		fprintf(logfile," done\n"); fflush(logfile);
	}*/
}
