

 ! MADEIRA

 
     
      PROGRAM teste equdif
      IMPLICIT real (a-h,o-z)

	 real hx,ky,lambda,ll,u,x,rox,freq, Qx, Eyin, Ey, Dp, Pot			 
	 real e2linhas,e1linha, eo !, er
	 complex :: er=(1.591,-0.033)!permissividade relativa(er=e1linha-j*e2linhas)


	parameter (Eyin=1.6E4)	 
	!valor de entrada do campo elétrico, como esse valor depende da potência, frequência, eo e e2linhas...
	! ele é calculado com base nesses valores na equação abaixo.		 
	parameter(freq=2.45E9) !frequência (Hz)
	parameter(rox=721)	!densidade (kg/m^3)		 art.numerical simulation...
	parameter(Cx=1126.8)	 !calor específico (J/kg*k^-1)
	parameter(CT=0.259)	  !condutividade térmica (W/mK)	(no art. é o K)
	parameter(eo=8.854E-12)  !permissividade no vácuo (F/m) 
	parameter(e1linha=1.591)	 !Constante dielétrica
	parameter(e2linhas=0.033)	 ! fator de perda
      parameter(xl=0.11)   ! comprimento máximo (m)
	parameter(T=20)  ! tempo máximo  (s)
	parameter(mx=100) !m do metodo
	parameter(Nt=100)  ! n do método
	parameter(v=2.999792457E8)	  !velocidade de fase  
	!parameter (Pot=1000)

	 

       dimension Ey(mx), ll(mx),u(mx),z(mx)
	 dimension Qx(mx), Fx(mx), w(mx),winic(mx)

	 open(6, file="saida.dat")

	   alfa=CT/(rox*Cx)
	  !difusividade térmica (m^2/s). Valor retirado do art. Numerical Simulation...p.16
	  write(*,*) 'alfa',alfa



	  hx=xl/mx
	  write(*,*) 'hx',hx

	  ky=T/Nt
	  write(*,*) 'ky',ky

	  lambda=((alfa)*ky)/(hx**2)

	  write(*,*) 'lambda',lambda



	   !penetração 

	   a=SQRT((e1linha*SQRT((1+(e2linhas/e1linha)**2))-1)/2)
                   
	   Dp=1/((2*3.1416*freq/v)*a)	   

	   write(*,*) 'Dp',Dp

	  
	  do i=1,mx-1,1				 
            
		 winic(i)=298.0   ! temp inicial

	      enddo

	 write(*,*) 'winic',winic

	!******** Campo elétrico de entrada. Equação original Pot=2*3.1416*freq*eo*e2linhas*|Eyin|**2


	  ! Eyin=SQRT(Pot/(2*3.1416*freq*eo*e2linhas))
	   !write(*,*) 'Eyin', Eyin


	  Pot=(2*3.1416*freq*eo*e2linhas)*Eyin**2
	   write(*,*) 'Pot', Pot


	   	  
		
	!************** iniciar loop do tempo de exposição


	    DO tj=1,T,1
	
     		 do i=1,mx-1,1           
		      w(i)=winic(i)   

	        enddo
     	    write(*,*) 'w',w
     
          do i=1, mx-1,1
		
		Ey(mx)=0.0

          enddo
	  

	    do i=1,mx-1,1
		  

	   Ey(i)=Eyin*sin((3.1416*i*hx)/xl)*sin(2*3.1416*freq*tj)	 
           !variação do campo elétrico

	   	  
	   enddo
	 write(*,*) 'Ey',Ey


	   do i=1,mx-1,1

	   Qx(i)=2*3.1416*freq*eo*er*(e2linhas/e1linha)*Ey(i)**2  
	   ! Q da equação do artigo
	  
          
	   enddo
	 !write(*,*) 'Qx',Qx


	    do i=1,mx-1,1
   
          Fx(i)=Qx(i)/(rox*Cx)  !Essa é a função depois da igualdade
          

	    enddo

		 write(*,*) 'Fx',Fx
  

	      


	    
	  

	 ll(1)=1+lambda
	 u(1)=(-lambda)/(2*ll(1))	  ! passo 3
	  

	  
	    do i=2,mx-2,1
	    
	       ll(i)=1+lambda+(lambda*u(i-1)/2)
	       									! passo 4
		   u(i)= -lambda/(2*ll(i))

	    enddo

        ll(mx-1)=1+lambda+(lambda*u(mx-2)/2) !passo 5 algoritmo pag 735 burden

	

	    
	    do j=1,Nt	   !passo 6
	          
			  	     
	z(1)=(((1-lambda)*w(1))+((lambda/2)*w(2))+ky*Fx(j))/(ll(1))	 			 
	 !passo 7

	          do i=2,mx-1


      z(i)=(((1-lambda)*w(i))+(lambda/2)*(w(i+1)+w(i-1)+z(i-1))+
     * ky*Fx(i))/(ll(i))
	!passo 8

!	write(*,*) 'z(i)',z



	          enddo
	     
	      w(mx-1)=z(mx-1)	 !passo 9


	          do i=mx-3,3,-1

	           
			  w(i)=z(i)-(u(i)*w(i+1))	 !passo 10
	          
			  enddo



	       do i=3,mx-3,1           
		       winic(i)=w(i)   ! temp inicial

	        enddo


	     
	          

	         do i=1,mx-1

	           x=i*hx

	         enddo

	    enddo !j

	ENDDO  !FECHAR TEMPO DE EXPOSIÇÃO



		 do i=3, mx-3,1
		 
		 write(6,*)i*hx, Ey(i), Qx(i), winic(i)

		 enddo
	 	 
		  
	 !ENDDO  !FECHAR TEMPO DE EXPOSIÇÃO	  
		  
		  	     
	    close(6)


        

   !************ fechaR LOOP DE TEMPO DE EXPOSIÇÃO

    	END PROGRAM !equação diferencial

