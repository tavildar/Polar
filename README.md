Polar codes 
===================

This repository provides C and MATLAB implementations for polar codes.

> For the seminal work on polar codes, please refer to: **Erdal Arikan**, "Channel polarization: A method for constructing capacity-achieving codes for symmetric binary-input memoryless channels",  http://arxiv.org/abs/0807.3917 

Overview of what is provided
----------

 - Encoding for polar codes

 - Decoding for polar codes including
	 - Successive cancellation (SC) decoding (See [Arikan])
	 - Successive cancellation list (SCL) decoding (See [Tal])
	 - LLR based SCL decoding (See [Stimming])
	 
 - Code construction: Bhattacharya parameter based construction 

 - Binary AWGN simulations

Decoding performance
------


Runtime performance C and MATLAB
-----
Runtime is mainly dominated by the decoder. The run time comparison for rate 1/2 code is as follows (run on a single macbook pro 2015):
<table>
<caption> Comparison of number of iterations per second </caption>
  <tr align="center">
    <th>Parameters </th>
    <th>PolarC </th>
    <th>PolarM</th>
	<th>Speedup C/M</th>
  </tr>
  <tr align="center">
    <td>N = 2048, L = 1</td>
    <td>250</td>
    <td>51 </td>
    <td>5x</td>
  </tr>
  <tr align="center">
    <td>N = 2048, L = 4</td>
    <td>67 </td>
    <td>3.2</td>
    <td>20x</td>
  </tr>
  <tr align="center">
    <td>N = 2048, L = 32</td>
    <td>10.1</td>
    <td>0.65</td>
    <td>15x</td>
  </tr>
  <tr align="center">
    <td>N = 512, L = 1</td>
    <td>890</td>
    <td>204</td>
    <td>4.5x</td>
  </tr>
  <tr align="center">
    <td>N = 512, L = 4</td>
    <td>287</td>
    <td>14</td>
    <td>20x</td>
  </tr>
  <tr align="center">
    <td>N = 512, L = 32</td>
    <td>44</td>
    <td>2.8</td>
    <td>16x</td>
  </tr>
</table>


Code Interface
------

The key code is in class PolarCode (C and MATLAB). The interface to this class is as follows:

 - **constructor**
	 - Input: N or n, K, epsilon, CRC 
	 - Output: None 
	 - Function: create polar code with specified parameters. 
		 - Bhattacharya parameter based construction is used.

 - **encode** (MATLAB code modified [Pfister])
	 - Inputs: info bits (length = K)
	 - Output: coded bits (length = N)

 - **decode_SC_P1** (only for MATLAB - modified from [Pfister]))
	 - Input: vector of probability of a bit being  '1' (length = N)
	 - Output: decoded info bits (length = K)
	 - Function: SC decoder  (see [Arikan])

 - **decode_SCL_P1**
	 - Input: p1, p0, L 
		 - p1 = vector of probability of output given bit = 1, 
		 - p0 = vector of probability of output given bit = 0,   
		 - L = list size
	 - Output: decoded info bits
	 - 	Function: SCL decoder (see [Tal])

 - **decode_SCL_LLR**
	 - Input: LLR, L
		 - LLR = vector of log likelihood ratios 
		 - L = list size
	 - Output: decoded info bits
	 - Function: LLR based SCL decoder (see [Stimming])

 - **get_bler_quick** (helper function)
	 - Input: EbNo_vec, list_size_vec
	 - Output: BLER estimate
	 - Function: get BLER estimate. 
		 
> Note that for speed up of simulation, get_bler_quick assumes that if a given run is decoded for a lower EbNo value, then it will be decoded for a higher EbNo value. 


References
---------

 1. **E. Arikan**, "Channel polarization: A method for constructing capacity-achieving codes for symmetric binary-input memoryless channels",  http://arxiv.org/abs/0807.3917 
 2.  **Ido Tal** and **Alexander Vardy**, 	"List Decoding of Polar Codes", https://arxiv.org/abs/1206.0050
 3.  **Alexios Balatsoukas-Stimming**, **Mani Bastani Parizi**, **Andreas Burg**, "LLR-based Successive Cancellation List Decoding of Polar Codes", https://arxiv.org/abs/1401.3753
 4. **Henry D. Pfister,** "A Brief Introduction to Polar Codes", http://pfister.ee.duke.edu/courses/ecen655/polar.pdf

