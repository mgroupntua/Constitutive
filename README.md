![alt text](http://mgroup.ntua.gr/wp-content/uploads/2018/05/MGroup52.png "MGroup")
# Constitutive

MSolve library that includes standard material models (phenomenological laws) as well as implementations of 
multiscale homogenization schemes for the modeling of microscopically heterogeneous materials.

## Features

### 1. Material library
The material library implemented in msolve consits of the following material laws:

1. Elastic material for 2D and 3D problems
2. Elastic material implementation suitable for shell finite element formulations
 based on the Kirchhoff-Love or the Reisner Mindlin hypotheses
3. Von-Mises elatoplastic material for 3D problems
4. Bond slip traction separation law for cohesive finite elements

### 2. Multiscale analysis

#### Main capabilities 
1. Concurrent (FE2) Multiscale analysis of microscopically heterogeneous structures
2. Stress strain history analysis of heterogenous material models

#### Two scale transitions implementations
1. small and finite strains implementations for 3D problems
2. small strains 3D to 2D homogenizations
3. small strains 3D to Kirchhoff Love shell
4. thermal problems

#### Model General properties
Classes of the Msolve Constitutive project are also used to assign general model properites:
1. for dynamic analysis
2. steady state thermal problems

## Installation instructions
You can choose either to clone the solution or downloads it as a zip file.

### Clone solution
1. Under the repository name, click **Clone or Download** option.

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/CloneOrDownload.png "1")

2. In the popup appearing choose the **Use HTTPS** option.

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/2.png "2")

3. Use the ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/3.png "3") to copy the link provided.

4. Open Visual Studio. In Team Explorer window appearing in your screen under Local Git Repositories click the **Clone** option. If Team Explorer window is not visible you can enable in View -> Team Explorer

  ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/4.png "4")
  
5. In the text box appearing paste the link.

 ![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/5.png "5")

6. Click clone and Visual Studio will automatically download and import **MGroup.Constitutive**


### Download as ZIP
1. Under the repository name, click **Clone or Download** option

![alt text](https://github.com/mgroupntua/MSolve.Edu/blob/master/Images/CloneOrDownload.png "1")

2. Click **Download ZIP** option. **MGroup.Constitutive** will be downloaded as a ZIP file.

3. Extract the ZIP file to the folder of choice.

4. Double click on **MGroup.Constitutive.sln** file to open the code with Visual Studio

  



