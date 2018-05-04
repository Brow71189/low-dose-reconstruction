#ifndef ZIPPEDIMG_HPP
#define ZIPPEDIMG_HPP

//fill zipImg with compressed data of qImg
long write_zipImg(	unsigned short* zipImg,	unsigned short* qImg);
//fill qImg with decompressed data from zipImg
void read_zipImg(	unsigned short* zipImg,	unsigned short* qImg);

//retrievs the value from the current zip_img_sfP and advances the pointer to the next subframe
//WARNING the zip_img_sfP pointer has to be reset manually
//there is also no check for the last valid subframe
//TODO maybe use an entrance table for the different sf, and only read until value is found, but mostly zeros anyways
unsigned short look_up_sf(int indpp);


#endif
