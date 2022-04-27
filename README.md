# SenK

SenK is a high-performance linear solver library in C++. The most important feature of SenK is its ILUB preconditioner. 

## ILUB precondtioner

$$
\mathbf{A} = \left[
	\begin{array}{cc|cc|cc}
	● &   &   &   & ● &    \\
	  & ● &   &   &   & ●  \\\hline
	  &   & ● &   &   &    \\
	● &   &   & ● & ● &    \\\hline
	  &   &   &   & ● &    \\
	  &   &   &   &   & ●  \\
	\end{array}
\right] \quad→\quad
\mathbf{L} = \left[
	\begin{array}{cc|cc|cc|cc}
	● &   &   &   &   &    \\
	● & ● &   &   &   &    \\\hline
	● & ● & ● &   &   &    \\
	● & ● & ● & ● &   &    \\\hline
	  &   &   &   & ● &    \\
	  &   &   &   & ● & ●  \\
	\end{array}
\right],
\mathbf{U} = \left[
	\begin{array}{cc|cc|cc|cc}
	● & ● &   &   & ● & ●  \\
	  & ● &   &   & ● & ●  \\\hline
	  &   & ● & ● & ● & ●  \\
	  &   &   & ● & ● & ●  \\\hline
	  &   &   &   & ● & ●  \\
	  &   &   &   &   & ●  \\
	\end{array}
\right]
$$

## Sample

A sample source code to execute ILUB precondtioned`sample.cpp` 

```c++
#include "senk.hpp"
int main(int argc, char **argv) {
  
}
```

