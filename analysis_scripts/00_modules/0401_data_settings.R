

## settings
if(TRUE) {
  
  project = 'Mel.20240206'
  
  output = project
  
  model.type = 'xgboost' 
  
  threads = 12 
  
  ## ---------------------------------------------------------------
  
  folder = '20240303a_archive_CV10'
  
  in.path = paste0('results/20240206/',folder)
  print(in.path)
  
  ## ---------------------------------------------------------------
  
  ## features 
  clinicvars = c()
  clinicvars = c(
    'Age',  
    'Sex', 
    'BMI',  
    'Melanoma.subtype.clp', 
    'NLR.recalc', 
    'eos.count.pre', 
    'Albumin.pre', 
    'LDH.pre',  
    'Mets.per.AJCC.01'
  )
  
  print(paste0('clinicvars = ', clinicvars))

  
  ## ---------------------------------------------------------------
  
  feat.orders = c('First_Order', 'Second_Order')
  ct.tissues = list(
    'AL' = 'Adrenal',
    'LN' = 'LN',
    'Liver_L' = 'Liver_L',
    'Lung_L' = 'Lung_L',
    # 'SL' = 'Splenic',
    'STL' = 'SoftTissue',
    'WCT' = 'Normalization'
  )
  
  # AL : Adrenal Lesion
  # Liver_L : Liver Lesion
  # LN : Lymph Node
  # Lung_L : :Lung Lesion
  # SL : Splenic Lesion
  # STL : Soft Tissue Lesion
  
  mri.tissues = list(
    'ED' = 'Edema',
    'EN' = 'Enhancement',
    'HEM' = 'Hemorrhage',
    'N' = 'Necrosis'
  )
  
  mri.types = list(
    'FLAIR' = 'FLAIR',
    'POST_T1' = 'POST_T1',
    'PRE_T1' = 'PRE_T1'
  )
  
  
  ## ---------------------------------------------------------------
  
  cohorts = c('all')
  runs = c(1)
  
  ## ---------------------------------------------------------------
  
  ## for caret...
  preproc.steps = c('scale', 'center','YeoJohnson')
  cor.abs.cutoff = 0.95 
  ratio = 0.8  
  cv.fold = 10 
  
  ## ---------------------------------------------------------------
  
  out.dir = file.path(my.cohort,'allsm')
  
  out.dir = file.path(model.type, out.dir)
  
  out.dir = paste0(out.dir,'/',clinicvars.set)
  
  if(!dir.exists(out.dir)) { dir.create(out.dir, recursive = T) }
  
  print(out.dir)
  
}