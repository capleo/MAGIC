# write control files for cambridge

#preliminary
options(stringsAsFactors = F)
library(qtl2)

#clean workspace, set wd
setwd("~/MAGIC/output/")

# write control files Joint
write_control_file("MAGIC.CAM.yaml",
                   crosstype = "genril9",
                   geno_file ="MAGIC_geno.csv",
                   founder_geno_file = "MAGIC_foundergeno9.csv",
                   gmap_file = "MAGIC_gmap_V4.csv",
                   pmap_file = "MAGIC_pmap_V4.csv",
                   pheno_file = "MAGIC_pheno.cam.Joint.all.csv",
                   geno_transposed = T, 
                   founder_geno_transposed = T,
                   crossinfo_file = "MAGIC_crossinfo.6.self.csv",
                   #geno_codes = genocodes,
                   #alleles = f_alls,
                   description = "MAGIC-CAM BLUPs joint model",
                   overwrite = T)

#write all temp
write_control_file("MAGIC.CAM.2122.yaml",
                   crosstype = "genril9",
                   geno_file ="MAGIC_geno.csv",
                   founder_geno_file = "MAGIC_foundergeno9.csv",
                   gmap_file = "MAGIC_gmap_V4.csv",
                   pmap_file = "MAGIC_pmap_V4.csv",
                   pheno_file = "MAGIC_pheno.cam.2122.all.csv",
                   geno_transposed = T, 
                   founder_geno_transposed = T,
                   crossinfo_file = "MAGIC_crossinfo.6.self.csv",
                   #geno_codes = genocodes,
                   #alleles = f_alls,
                   description = "MAGIC-CAM all BLUPs from 2021 and 2022",
                   overwrite = T)
