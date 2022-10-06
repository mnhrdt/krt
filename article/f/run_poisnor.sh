cat gbarbara.png                    |
	plambda 'x,gf dup vnorm /'  |  # compute normalized gradient (fwd diff)
	plambda 'x,db'              |  # compute divergence (backward diff)
	fft                         |  # fft
	plambda ':L /'              |  # anti-laplacian filter
	ifft                        |  # ifft
	blur C 1                    |  # remove very low frequencies
	qauto -p 0                  |  # color balance
	cpu
