import os, sys
import glob
import re
import astropy as a
import ccdproc
import numpy as n
import time

#PARAMETERS
Telescope = "LOT"
Camera = ""
object_name = "NGC7006"
DARK_DURATION = 10

def classfile():
    return ["*.fts", "*.fit", "*.fits", "*.tif", "*.tiff", "*.NEF"]
cf = classfile()
date = [time.localtime()[t] for t in range(0, len(time.localtime())-3)]
#print(date)

def class_img():
    return ["BIAS", "DARK", "FLAT"]
    
def class_filter():
    return ["gp", "rp", "ip"]

def create_direc_BDF():
    direc = os.listdir()
    print(direc)
    for ci in class_img():
        if ci in direc:
            print(f"Directory {ci} exists!")
        else:
            print(f"Directory {ci} doesn't exist! We will create {ci}")
            time.sleep(1)
            os.mkdir(ci)
            print(f"Directory {ci} is already set up!")   

raw_path = os.path.realpath(sys.argv[0])
raw_direc = os.path.dirname(raw_path)
#print(raw_direc) #Checkpoint

raw_file_list = glob.glob("./*", recursive=True)
#print(raw_file_list) #Checkpoint

def moving_BDF():
    for st in cf:
        root_direc = os.path.join(raw_direc, "**", st)
        #print(root_direc) #Checkpoint

        root_file_list = glob.glob(root_direc, recursive=True)
        #print(root_file_list) #checkpoint
        for ci in class_img():
            count, mcount = 0, 0
            print(f"Searching {ci}.{st}*:")
            for rfl in root_file_list:
                imgtype = rfl.split("/")
                if ci.lower() in imgtype[-1].lower():
                     print(imgtype[-1])
                     if os.path.dirname(rfl) == raw_direc + f"/{ci}":
                         print(f"File {imgtype[-1]} is already at directory {ci}")
                     else:
                         print(f"File {imgtype[-1]} is not at directory {ci}")
                         time.sleep(0.1)
                         print(f"Move {imgtype[-1]} to directory {ci}")
                         os.replace(rfl, raw_direc + f"/{ci}" + f"/{imgtype[-1]}")
                         time.sleep(0.1)
                         mcount += 1
                     #print(os.path.dirname(rfl))
                     #print(raw_direc + f"/{ci}")
                     count += 1
            print("===found {} files===\n===move {} files to workspace===\n".format(count, mcount))    
            time.sleep(0.1)
            
def search_header():
    while True:
    	inc = 0
    	ftls = []
    	for st in cf:
    	    root_direc = os.path.join(raw_direc, "**", st)
    	    root_file_list = glob.glob(root_direc, recursive=True)
    	    for rfl in root_file_list:
    	        print(f"[{inc}]", rfl)
    	        inc += 1
    	        ftls.append(rfl)
    	nl = [str(n) for n in range(0, len(ftls)+2)]
    	ccinc = input("Which header do you want to see? ([n]) Ans: ")
    	if ccinc not in nl:
    	    print("Invalid keyword! Try again!")
    	    continue
    	elif int(ccinc) < 0 or int(ccinc) >= len(ftls):
    	    print("Outer Index! Please try again!")
    	    continue
    	else:
    	    fts = a.io.fits.open(ftls[int(ccinc)])
    	    ftsh = fts[0].header
    	    fts.close()
    	    #print(ftsh)
    	    for kf in ftsh:
    	        print(kf, ftsh[kf])
    	    print("IMAGETYP header value:", ftsh["IMAGETYP"])
    	break

def stacking_BDF():
    for ci in class_img():
        root_direc = raw_direc + f"/{ci}"
        #print(root_direc) #Checkpoint
        print(f"Prepare to stack {ci}")
        time.sleep(1)
        ccdp_img_col = ccdproc.ImageFileCollection(root_direc)
        print(ccdp_img_col)
        
        #Stacking BIAS
        if ci == class_img()[0]:
            select_bias = ccdp_img_col.files_filtered(
              IMAGETYP= ci,
              TELESCOP="LOT",
              CAMERA="Andor DZ936 S/N",
              XBINNING=1,
              GAIN="2", #be careful of the data type LOL, our fits use string her
              include_path=True)
            for bias in select_bias:
                print(bias)
            print(f"Combine {ci}...")
            time.sleep(1)
            combine_bias = ccdproc.combine(select_bias,
              unit="adu",
              method='average',
              sigma_clip=True,
              sigma_clip_low_thresh=3,
              sigma_clip_high_thresh=3,
              sigma_clip_func=n.ma.median,
              sigma_clip_dev_func=a.stats.mad_std,
              mem_limit=350e6,
              dtype="float32")
            print(f"=== Finish combining {ci} ===")
            time.sleep(1)
            print(f"Store the {ci} combined file...")
            time.sleep(1)
            combined_bias_dir = os.path.join(root_direc, f"master_{ci}.fits")
            combine_bias.write(combined_bias_dir, overwrite=True)
            print(f"=== Finish storing {ci} combined file! ===")
            time.sleep(1)
        
        #Stacking DARK    
        elif ci == class_img()[1]:
            print(f"Find the combined {class_img()[0]} file...")
            BIAS_direc = os.path.join(raw_direc + f"/{class_img()[0]}", "**", classfile()[2])
            BIAS_combined_file = glob.glob(BIAS_direc, recursive=True)
            BIAS_combined_read = a.nddata.CCDData.read(BIAS_combined_file[0], unit = "adu")
            #print(BIAS_combined_read) #Checkpoint
                 
            select_dark = ccdp_img_col.files_filtered(
            IMAGETYP= ci,
            TELESCOP="LOT",
            CAMERA="Andor DZ936 S/N",
            XBINNING=1,
            GAIN="2",
            EXPTIME=DARK_DURATION, #be careful of the data type LOL, our fits use string here
            include_path=True)
            time.sleep(1)
            
            print(f"Finding {ci} LIST")
            print(f"======{ci} LIST======")
            for dark in select_dark:
                print(dark)
            print("======LIST END HERE======")
            time.sleep(1)
            
            #Workspace for subtracted BIAS
            print("Searching \"subbias\" directory")
            time.sleep(1)
            subdirec = os.listdir()
            #print(subdirec) #Checkpoint
            if "subbias" not in subdirec:
                print("There isn't \"subbias\" directory in our workspace")
                time.sleep(0.5)
                print("Creating \"sub\" directory for sub tracted files for the first time")
                os.mkdir("subbias")
                time.sleep(1)
                print("Already created the \"subbias\" directory!")
                time.sleep(1)
            else:
                print("There is \"subbias\" directory in our workspace")
            print("===Finish setting===")    
            time.sleep(1)
            
            #Subtract BIAS
            print(f"Subtracting {class_img()[0]} from {ci} images...")
            time.sleep(1)
            for dark in select_dark:
                sbf = "subbias_" + dark.split("/")[-1]
                DARK_read = a.nddata.CCDData.read(dark, unit = "adu")
                DARK_calibrated = ccdproc.subtract_bias(DARK_read, BIAS_combined_read)
                DARK_calibrated.data = DARK_calibrated.data.astype("float32")
                print(f"Store \"{sbf}\"...")
                subbias_file = os.path.join(raw_direc + "/subbias", sbf)
                DARK_calibrated.write(subbias_file, overwrite = True)
            print("===Finish Subtracting===")
            time.sleep(1)
            
            print(f"Combine {ci}...")
            time.sleep(1)
            combine_dark = ccdproc.combine(select_dark,
              unit="adu", 
              method='average',
              sigma_clip=True, 
              sigma_clip_low_thresh=3, 
              sigma_clip_high_thresh=3, 
              sigma_clip_func=n.ma.median, 
              sigma_clip_dev_func=a.stats.mad_std, mem_limit=350e6, dtype="float32")
            print(f"=== Finish combining {ci} ===")
            time.sleep(1)
            print(f"Store the {ci} combined file...")
            time.sleep(1)
            combined_dark_dir = os.path.join(root_direc, "master_{}_{}-s.fits".format(ci,
              str(DARK_DURATION)))
            combine_dark.write(combined_dark_dir, overwrite=True)
            print(f"=== Finish storing {ci} combined file! ===")
            time.sleep(1)
        
        #Stacking FLAT
        elif ci == class_img()[2]:
            print(f"Find the combined {class_img()[0]} and {class_img()[1]} file...")
            BIAS_direc = os.path.join(raw_direc + f"/{class_img()[0]}", "**", classfile()[2])
            BIAS_combined_file = glob.glob(BIAS_direc, recursive=True)
            BIAS_combined_read = a.nddata.CCDData.read(BIAS_combined_file[0], unit = "adu")
            DARK_direc = os.path.join(raw_direc + f"/{class_img()[1]}", "**", classfile()[2])
            DARK_combined_file = glob.glob(DARK_direc, recursive=True)
            DARK_combined_read = a.nddata.CCDData.read(DARK_combined_file[0], unit = "adu")
            time.sleep(1)
            
            for filt in class_filter():
                print(f"===calibrate filter: {filt}===")
                select_flat = ccdp_img_col.files_filtered(
                  IMAGETYP= ci,
                  TELESCOP="LOT",
                  CAMERA="Andor DZ936 S/N",
                  XBINNING=1,
                  GAIN="2", #be careful of the data type LOL, our LOT fits use string here. e/adu
                  FILTER=f"{filt}_Astrodon_2019",
                  include_path=True)
            
                print(f"Finding {ci}:{filt} LIST")
                print(f"======{ci}:{filt} LIST======")
                for flat in select_flat:
                    print(flat)
                print("======LIST END HERE======")
                time.sleep(1)
             
                #Subtract BIAS
                print(f"Subtracting {class_img()[0]} from {ci}:{filt} images...")
                time.sleep(1)
                for flat in select_flat:
                    sbf = "subbias_" + flat.split("/")[-1]
                    FLAT_read = a.nddata.CCDData.read(dark, unit = "adu")
                    FLAT_subbias = ccdproc.subtract_bias(FLAT_read, BIAS_combined_read)
                    FLAT_calibrated = ccdproc.subtract_dark(FLAT_subbias, DARK_combined_read,
                      exposure_time = 'exptime', 
                      exposure_unit = a.units.second, 
                      scale = True)
                    FLAT_calibrated.data = FLAT_calibrated.data.astype("float32")
                    print(f"Store \"{sbf}\"...")
                    subbias_file = os.path.join(raw_direc + "/subbias", sbf)
                    FLAT_calibrated.write(subbias_file, overwrite=True)
                print("===Finish Subtracting===")
                time.sleep(1)
                
                print(f"Combine {ci}...")
                def inv_median(a):
                    return 20000/n.median(a)
                time.sleep(1)
                combine_flat = ccdproc.combine(select_flat, 
                  unit="adu", 
                  method='median', 
                  scale=inv_median,
                  sigma_clip=True, 
                  sigma_clip_low_thresh=3, 
                  sigma_clip_high_thresh=3,
                  sigma_clip_func=n.ma.median, 
                  sigma_clip_dev_func=a.stats.mad_std,
                  mem_limit=350e6, 
                  dtype="float32")
                print(f"=== Finish combining {ci}:{filt} ===")
                time.sleep(1)
                print(f"Store the {ci}:{filt} combined file...")
                time.sleep(1)
                combined_flat_dir = os.path.join(root_direc, "master_{}_{}-f.fits".format(ci, filt))
                combine_flat.write(combined_flat_dir, overwrite=True)
                print(f"=== Finish storing {ci}:{filt} combined file! ===")
                time.sleep(1)       

def calibrate_img():
     #Find master file
     master_img = os.path.join(raw_direc, "**" ,"master" + classfile()[2])
     master_file = glob.glob(master_img, recursive = True)
     #Find image file
     image_img = os.path.join(raw_direc + "/Object" + f"/{object_name}", "**", classfile()[0])
     image_file = glob.glob(image_img, recursive = True)

     l_master = []
     cal_sciimg = []
     
     print(f"There are {len(image_file)} science files prepared to subtract")
     time.sleep(1)
     print(f"=== science files for target: {object_name} ===")
     for f_img in image_file:
          print(f_img.split("/")[-1])
     print("===END HERE===\n")
     time.sleep(1)
     
     #Sort BIAS, DARK, FLAT
     for ci in class_img():
         for mf_r in master_file:
             if ci.lower() in mf_r.lower():
                 l_master.append(mf_r)
             else:
                 continue
     
     print("MASTER FILES LIST")
     time.sleep(1)           
     print("=== master files ===") 
     for masimg in l_master:
         print(masimg.split("/")[-1])
     print("===END HERE===\n")
     time.sleep(1)
     
     for sci in image_file:
         sci_subbias, sci_subdark, sci_subflat = [], [], []
         ccd_sci = a.nddata.CCDData.read(sci, unit = "adu")
         
         for mf in l_master:
             if class_img()[0].lower() in mf.lower():
                 print("===subtract {} within sub{}===".format(sci.split("/")[-1], class_img()[0].lower()))
                 ccd_mbias = a.nddata.CCDData.read(mf, unit = "adu")
                 subbias = ccdproc.subtract_bias(ccd_sci, ccd_mbias)
                 sci_subbias.append(subbias)
                 print("\n")
                 
             elif class_img()[1].lower() in mf.lower():
                 print("===subtract {} within sub{}===".format(sci.split("/")[-1], class_img()[1].lower()))
                 print(f"===subtract sub{class_img()[1].lower()}===")
                 ccd_mdark = a.nddata.CCDData.read(mf, unit = "adu")
                 subdark = ccdproc.subtract_dark(sci_subbias[0], ccd_mdark,
                             exposure_time = DARK_DURATION,
                             exposure_unit = a.units.second,
                             scale = True)
                 sci_subdark.append(subdark)
                 print("\n")
             
             elif class_img()[2].lower() in mf.lower():
                 print("===subtract {} within flat:{}===".format(sci.split("/")[-1], mf.split("/")[-1].replace("master_", "")))
                 mfi = mf.split("/")[-1]
                 sfi = sci.split("/")[-1]
                 for filt in class_filter():
                     if filt in mfi and filt in sfi:
                         print(f"Found masterfile: {mfi}!")
                         ccd_mflat = a.nddata.CCDData.read(mf, unit = "adu")
                         subflat = ccdproc.flat_correct(sci_subdark[0], ccd_mflat)
                         sci_subflat.append(subflat)
                         sci_subflat[0].data = sci_subflat[0].data.astype("float32")
                         
                         calibrated_science_dir = os.path.join(raw_direc + "/subbias", "sub{}".format(sci.split("/")[-1].replace(classfile()[0],classfile()[2])))
                         print(f"Store {calibrated_science_dir}")
                         sci_subflat[0].write(calibrated_science_dir, overwrite = True)
                         cal_sciimg.append(calibrated_science_dir)
                         print(f"Finish subtracting {sci}")
                         break
                     else:
                         continue
                         
                 print("\n")
                 
     time.sleep(1)
     print("Improving...")
     time.sleep(1)
     print("===calibrated scifile list===")
     for csi in cal_sciimg:
         print(csi.split("/")[-1])
     print("===END HERE===")
     time.sleep(1)
     
     #improve cal.sci_images
     for csi in cal_sciimg:
         ccd_csi = a.nddata.CCDData.read(csi, unit = "adu")
         ccd_csi.data = n.clip(ccd_csi.data, 1, 2 ** 16 - 1)
         ccd_csi.data[1] = 2 ** 16 - 1
         
         #Convert to LAcosmic sys.
         ccd_csi_elec = ccdproc.gain_correct(ccd = ccd_csi, gain = 2*a.units.electron/a.units.adu)
         #Remove cosmic ray with LAcosmic sys.
         ccd_csi_electron = ccdproc.cosmicray_lacosmic(
                              ccd=ccd_csi_elec,
                              satlevel=65535,
                              readnoise=8.5,
                              sigclip=7,
                              cleantype="medmask",
                              niter=4,
                              verbose=True)
         ccd_csi.data = ccd_csi_electron.data / 2
         ccd_csi.meta["BUNIT"] = "adu"
         ccd_csi.unit = "adu"
         
         #Bad pixel (column) map for the camera
         BPM = [{"x_pixel_coord_small": 951,#first bad column area
                 "x_pixel_coord_large": 953,
                 "y_pixel_coord_small": 1159,
                 "y_pixel_coord_large": 2047,},
                {"x_pixel_coord_small": 1565,#second bad column area
                 "x_pixel_coord_large": 1568,
                 "y_pixel_coord_small": 2013,
                 "y_pixel_coord_large": 2047,},
                {"x_pixel_coord_small": 1293,#second bad column area
                 "x_pixel_coord_large": 1295,
                 "y_pixel_coord_small": 1140,
                 "y_pixel_coord_large": 1149,},
                {}] #fourth bad column area and so on
         
         #Remove BPM
         ccd_csi_img_median_v = n.nanmedian(ccd_csi.data)
         bcm_pxrange = 2
         #Generate mask
         mask = n.ma.make_mask(ccd_csi.data,
                  copy=True,
                  shrink=True,
                  dtype=bool)
         mask[:,:]=False
         for bcc in BPM:
             if bcc != {} and len(bcc.keys()) == 4:
                 mask[int(bcc["y_pixel_coord_small"])
                       :int(bcc["y_pixel_coord_large"] + 1),
                        int(bcc["x_pixel_coord_small"])
                       :int(bcc["x_pixel_coord_large"] + 1)] = True
         
         masked_data = n.ma.masked_array(ccd_csi.data, mask = mask)
         csi_bpmf = ccd_csi.data.copy()

         for j in range(0,masked_data.shape[0]):
             for i in range(0,masked_data.shape[1]):
                 if mask[j,i] == True:
                    x1 = i - bcm_pxrange
                    x2 = i + bcm_pxrange + 1
                    y1 = j - bcm_pxrange
                    y2 = j + bcm_pxrange +1

                    if x1 < 0:
                        x1 = 0
                    if x2 > masked_data.shape[1]:
                        x2 = masked_data.shape[1]
                    if y1 < 0:
                        y1 = 0
                    if y2 > masked_data.shape[0]:
                        y2 = masked_data.shape[0]

                    replaced_data = n.ma.median(masked_data[y1:y2,x1:x2])
            
                    if n.isnan(float(replaced_data)) == True:
                        replaced_data = ccd_csi_img_median_v
            
                    csi_bpmf[j,i] = replaced_data
         
         ccd_csi.data = csi_bpmf.data
         
         #Create direc for saving improved_csi
         direc = os.listdir()
         if "improv" not in direc:
             print("create directory:", "improv")
             time.sleep(1)
             os.mkdir("improv")
             print("directory [improv] is created!")
             time.sleep(1)
         print("Store {}".format(csi.split("/")[-1].replace("sub","")))
         time.sleep(1)
         improved_ccd_csi = os.path.join(raw_direc + "/improv", "improv_{}".format(csi.split("/")[-1].replace("sub",""))) 
         ccd_csi.write(improved_ccd_csi, overwrite = True)
         print("Finish storing {}".format(csi.split("/")[-1].replace("sub","")))
         time.sleep(1)
         
     print("Remove subbias...")
     time.sleep(1)
     os.chdir(raw_direc + "/subbias")
     subf = os.listdir(raw_direc + "/subbias")
     print(subf) 
     for dsf in subf:
         print("remove {}...".format(dsf.split("/")[-1]))
         os.remove(str(dsf.split("/")[-1]))
     os.chdir(raw_direc)
     os.rmdir("subbias")
     print("Already Removed!")
     
     time.sleep(1)
               
### Process End Here So Far###
#Combiner
#def combiner_for_sci():
    
       

