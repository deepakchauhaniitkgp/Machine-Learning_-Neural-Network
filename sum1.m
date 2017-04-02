function summ = sum1(a)
 length1 = size(a);
    length = length1(2);
    summ = 0;
    for i=1:length
    summ  = summ + a(i)*a(i);
    end
   
end