function return_data = H_exact(params)
%H_EXACT Exact calculation of the convolution kernel as derived at the IPAM
%hackathon.
%
%last updated 03/10/26 by Adam Petrucci

    y = linspace(-pi,pi,params.N_y);
    inp = periodic_influence(0,y,2*pi);
    return_data = [h11(inp,params.del),zeros(size(inp)),...
                   zeros(size(inp)),h22(inp,params.del)];

end

function return_data = h11(r,del)
%H11 The (1,1) element of the H matrix
%
%last updated 03/10/26 by Adam Petrucci
    
    L = 2*pi;
    num = 4*L;
    denom = L^2+4*jap_brac(r,del).^2;
    return_data = (num./denom);

end

function return_data = h22(r,del)
%H22 The (2,2) element of the H matrix
%
%last updated 03/10/26 by Adam Petrucci

    L = 2*pi;
    jb = jap_brac(r,del);

    t1f1 = (del^2-r.^2)./(jb.^2) - 1;
    t1f2 = (2*L./(L^2+4*jb.^2));
    t2f1 = 1+((del^2-r.^2)./(jb.^2));
    t2f2 = atan(L./(2*jb))./jb;
    return_data = t1f1.*t1f2 + t2f1.*t2f2;

end