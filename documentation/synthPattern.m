function frame=synthPattern(x, amp, startIndex, width, spacing, gap)
% Creates synthetic pattern
%
% frame=synthPattern(x, amp, startIndex, width, spacing)

    frame= abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*2])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*3])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*4])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*5])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*6])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*7])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*8])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*9])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*10])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*11+gap])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*12+gap])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*13+gap])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*14+gap])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*15+gap])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*16+gap])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*17+gap])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*18+gap])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*19+gap])+...
    abs(rand(1)/100+1)*amp*gaussmf(x,[width,startIndex+spacing*20+gap]);

end