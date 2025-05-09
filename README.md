# **A NOVEL APPROACH FOR ASSESSING PILE BASE RESPONSE**
## **About**
This paper proposes a numerical approach to evaluate the base load-settlement response of cast-in-situ piles in sandy soils. 
The method assumes that the load response of the pile base can be modelled as the response of a rigid disk on a non-linear elastic half space following a hyperbolic stress-strain model captured through a stepwise-linear incremental procedure.  
Furthermore, empirical expressions for the bearing capacity factor are derived for specific relative settlements, along with an approximate closed-form equation for the estimation of the load-settlement curve, representing easy-to-implement tools in practical pile design. 
The method performance is validated against calibration chamber tests, Finite Element Analysis and full-scale pile load tests.
## **Authors**
Raffaele Cesaro*, Rodrigo Salgado**, Monica Prezzi**, Raffaele Di Laora*, Alessandro Mandolini*

*_Department of Engineering, Università della Campania ‘Luigi Vanvitelli’, 81031 Aversa (CE), Italy_

*_Lyles School of Civil Engineering, Purdue University, West Lafayette, IN 47907-1284, USA_

email: raffaele.cesaro@unicampania.it, rodrigo@ecn.purdue.edu, mprezzi@ecn.purdue.edu, raffaele.dilaora@unicampania.it, alessandro.mandolini@unicampania.it 
## **Provided files**
The proposed numerical code is implemented in MATLAB. Here, the code is provided in two versions:
- _Numerical_model_Code.m_ : This file is designed to be used independently; the user can modify the code and change the input data directly within the MATLAB script.
  Running this file the load-settlement curve is obtained, along with contour plots that show the evolution of various parameters within the soil volume, corresponding to the last load level achieved during the analysis.
- _Numerical_model_ExcelCode.m_ & _Numerical model_Excel app.xlsm_ : These two files must be downloaded together and placed in the same folder.
  This version is designed for users who do not prefer to use MATLAB and are only interested in estimating the load-settlement curve. Therefore, you should open the Excel file, specify the MATLAB script path in cell A1 (remember to write the path without quotes),
  and enter the input data. Then, by clicking the 'Run the code' button, the MATLAB script will be launched automatically. The calculation will be completed in a few seconds and will be deemed finished when the message 'The specified MATLAB script has been launched' appears
  . Once the analysis is complete, the MATLAB code will have created a file called 'results.txt' in the same folder. The path to this file must be specified in cell A2 of the Excel file. By clicking the 'Import results' button, the results will be automatically
   imported into Excel, and the load-settlement curve will be updated. The Excel file also includes empirical correlations for simplified evaluation of the load-settlement curve and for estimating the mobilized resistance at a relative settlement _w/d_ of 5% and 10%.
  Note that when using these empirical correlations, it is not necessary to specify the maximum and minimum void ratios (_e<sub>max</sub>_ and _e<sub>min</sub>_).
  
  Remember that whenever the input data in the Excel file are changed, the file must be saved before running the analysis.

  N.B.: After downloading the Excel file, remember to enable macros and to unblock the file, right-click on the file in your downloads folder, select "Properties", and then check the "Unblock" box in the General tab. Click "Apply" and "OK" to save the changes. This is needed in case automatic restrictions are applied to downloaded files on your computer.

