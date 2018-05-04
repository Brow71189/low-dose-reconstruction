These are the ImageJ plugins required for the low-dose sampling and reconstruction process. To install them, copy them into the "plugins" folder of ImageJ and restart ImageJ. Then, on the first start of one of them, select "Compile and Run..." from the plugins menu in ImageJ and select the plugin you want to run. In a subsequent launch of ImageJ they should appear in the list of plug-ins.

Note that some of these plugins depend on the apache commons math module, which is also included in this repository as a separate folder (commons-math3-3.5). Make sure to copy the file "commons-math3-3.5.jar" from that folder to your ImageJ plugins folder in order to ensure everything works properly.

The plug-ins included here are:
    * Hex_Sampler
    * Hex_Collector
    * Hex_Magic
    * Dragon_Fly
    * Frame_Aligner
