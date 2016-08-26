function f=diff_bn_structs(satrec0,satrec)
f=[0,0];

         f(1)=abs(satrec.error-satrec0.error);
        f(2)=abs(satrec.satnum-satrec0.satnum);
       f(3)=abs(satrec.epochyr-satrec0.epochyr);
     f(4)=abs(satrec.epochdays-satrec0.epochdays);
         f(5)=abs(satrec.ndot-satrec0.ndot);
         f(6)=abs(satrec.nddot-satrec0.nddot);
         f(7)=abs(satrec.bstar-satrec0.bstar);
         f(8)=abs(satrec.inclo-satrec0.inclo);
         f(9)=abs(satrec.nodeo-satrec0.nodeo);
          f(10)=abs(satrec.ecco-satrec0.ecco);
         f(11)=abs(satrec.argpo-satrec0.argpo);
            f(12)=abs(satrec.mo-satrec0.mo);
            f(13)=abs(satrec.no-satrec0.no);
             f(14)=abs(satrec.a-satrec0.a);
          f(15)=abs(satrec.alta-satrec0.alta);
          f(16)=abs(satrec.altp-satrec0.altp);
    f(17)=abs(satrec.jdsatepoch-satrec0.jdsatepoch);

        
         
end
