Prerequisites for Deployment 

. Verify the MATLAB runtime is installed and ensure you    
  have installed version 8.5 (R2015a) or version 77.   

If there have no MATLAB runtime in your system, then please proceed according following instruction

Procedure to install:

1) Run the MCR_R2015a_win64_installer or MCR_R2015a_win32_installer according your system. 

   Download the required version of the MATLAB runtime for R2015a from the MathWorks Web site by navigating to

         URL: http://www.mathworks.com/products/compiler/mcr/index.html  

   For more information about the MATLAB runtime and the MATLAB runtime installer, see Package and Distribute in the MATLAB Compiler documentation in the MathWorks Documentation Center.    

   NOTE: You will need administrator rights to run MCRInstaller.

2) Put 'mclmcrrt77.dll' file into your System directory (By default, this is C:\Windows\System (Windows 95/98/Me), 

   C:\WINNT\System32 (Windows NT/2000), or C:\Windows\System32 (Windows XP, Vista, 7)). then reboot your computer.

3) Now you can install MyAppInstaller_web application. 

   Then run 'selectDistribution' file to get the project.(if still there have problem like, system need DllRegisterServer, then dowonload it and install it into your System)

Thank you so much.