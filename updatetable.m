function matrix = updatetable(s,a,v,h,angle,hh,totalT,matrix,M)
     matrix(1+s,1) = a;
     matrix(1+s,2) = v;
     matrix(1+s,3) = h;
     matrix(1+s,4) = angle;
     matrix(1+s,5) = hh;
     matrix(1+s,6) = totalT;
     matrix(1+s,7) = M;
end
