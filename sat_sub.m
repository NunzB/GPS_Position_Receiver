function[sub_const_matrix] = sat_sub(dim_sub_const,sat_matrix,id_min,r)

        sub_const_matrix = zeros(dim_sub_const,3);
    
        for i=1:dim_sub_const
        for k=1:r
        if sat_matrix(k,4) == id_min(i)
           sub_const_matrix(i,1) = sat_matrix(k,1) ;
           sub_const_matrix(i,2) = sat_matrix(k,2);
           sub_const_matrix(i,3) = sat_matrix(k,3);    
        end
        end
        end
         
end
