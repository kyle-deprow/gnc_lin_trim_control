function [ derivative ] = first_order_derivative( array, h )

l = size(array,2);
derivative = zeros(l,1);
derivative(1) = (-3*array(1)+4*array(2)-array(3))/(2*h);
for i = 2:l-1
    derivative(i) = (array(i+1)-array(i-1))/(2*h);
end
derivative(l) = (array(l-2)-4*array(l-1)+3*array(l))/(2*h);

end

