# Applying-dimensionality-reducion-techniques-for-measured-caridac-blood-flow-
# About our Dataset :
* we have patient-specific blood flow information measured using four-dimensional phase-contrast magnetic resonance imaging, also known as 4D-PC-MRI i.e a 3D blood flow information.
* 197 attributes and 122 rows, each row represents characteristics of a patient 
* Attributes can be divided into 5 types 
  * Pressure related attributes
  * Velocity related attributes
  * Diameter related attributes 
  * Vortex related attributes
  * Flow jet related attributes 
# Task :
* Cluster patients into 3 cohorts 
  * Healthy patients 
  * Patients suffering from Bicuspid aortic valve disease 
  * Tetralogy of fallot.
# Methods used :
1.	Principal component analysis 
2.	T-distributed stochastic neighboring embedding 
3.	Uniform manifold approximation and projection 
4.	Multi dimensional scaling 
# Evaluation :
* Qualitative Evaluation :
  * Silhouette coefficient 
  * Correlation coefficient 
  * Connectivity validity measure
* Quantitative:
   * User study on 40 samples

additionally we have implemented our pre trained parameters (parameters of dimensionality techniques used on blood flow dataset) on a mice-protietn dataset, to check the cohort seperability

