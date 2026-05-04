## FFTW3_InterfaceNd

The SCC::FFTW3_InterfaceNd classes N = 1,2,3 are classes that facilitate the use of the fast discrete Fourier transforms provided by FFTW3. The interfaces provide scalings and coefficient organization that is useful when implementing spectral methods, e.g. methods that use a basis of discrete Fourier modes of both positive and negative wave numbers.

The MKL version of the FFTW3 library does not implement the real-to-real (r2r) transforms required by SCC_fftw3_sin2d.h and SCC_fftw3_sin3d.h. Programs will successfully link to the MKL FFTW library, but when  routines are used an exception will be thrown with a message similar to 


    XXX Error : required fftw_r2r_2d function not available
    XXX in FFTW library used (likely MKL)



### Prerequisites
C++11, FFTW3 installed, SCC::DoubleVectorNd, SCC::GridFunctionNd
### Versioning
Release : 1.0.10
### Date
Aug. 23, 2021
### Authors
Chris Anderson
### License
GPLv3  For a copy of the GNU General Public License see <http://www.gnu.org/licenses/>.
### Acknowledgements









