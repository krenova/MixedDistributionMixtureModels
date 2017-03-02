
# 1age: age in years
# 2sex: sex (1 = male; 0 = female)
# 3cp: chest pain type
        # -- Value 1: typical angina
        # -- Value 2: atypical angina
        # -- Value 3: non-anginal pain
        # -- Value 4: asymptomatic
# 4trestbps: resting blood pressure (in mm Hg on admission to the hospital)
# 5chol: serum cholestoral in mg/dl
# 6fbs: (fasting blood sugar > 120 mg/dl)  (1 = true; 0 = false)
# 7restecg: resting electrocardiographic results
        # -- Value 0: normal
        # -- Value 1: having ST-T wave abnormality (T wave inversions and/or ST 
                    # elevation or depression of > 0.05 mV)
        # -- Value 2: showing probable or definite left ventricular hypertrophy
                    # by Estes' criteria'
# 8thalach: maximum heart rate achieved
# 9exang: exercise induced angina (1 = yes; 0 = no)
# 10oldpeak = ST depression induced by exercise relative to rest
# 11slope: the slope of the peak exercise ST segment
        # -- Value 1: upsloping
        # -- Value 2: flat
        # -- Value 3: downsloping
# 12ca: number of major vessels (0-3) colored by flourosopy
# 13thal: 3 = normal; 6 = fixed defect; 7 = reversable defect
# 14num: diagnosis of heart disease (angiographic disease status)
        # -- Value 0: < 50% diameter narrowing
        # -- Value 1: > 50% diameter narrowing
        # (in any major vessel: attributes 59 through 68 are vessels)
		

output <- MFMM(cleveland_sub[,-14], K = 2, est_cov = TRUE, attempts = 10)

#accuracy
sum(apply(output$data_prob,1,function(x)which(x==max(x)))==((cleveland_sub$num>0)*1+1))

#confusion matrix
table(apply(output$data_prob,1,function(x)which(x==max(x))),((cleveland_sub$num>0)*1+1))