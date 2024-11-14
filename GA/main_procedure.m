function [record, tbest] = ...
    main_procedure(Pc, Pm, Pn, max_iter, pop_size, power, kfzd, t, lower, upper, pr, ip)
if kfzd(ip) == 1
    dim = t(ip);
    chromo_size = zeros(dim, 1);
    for ii = 1:1:t(ip)
        chromo_size(ii) = ceil(log2((upper(ip) - lower(ip) + 1) / 0.1 + 1));
    end
    pop = cell(1, dim);
    pop_int = zeros(pop_size, dim);
    for iii = 1:dim
        pop{1, iii} = initilize_pop(pop_size, chromo_size(iii));
    end
    record = zeros(1, max_iter);
    best = inf;
    for iter = 1:1:max_iter
        for il = 1:1:dim
            pop_int(:,il) = ...
                round(pop_decode(pop{1,il}, pop_size, chromo_size(il), upper(ip), lower(ip)));
        end
        Fit = zeros(pop_size, 1);
        Fit_ck = zeros(pop_size, 1);
        for io = 1:1:pop_size
            T = zeros(1,96);
            index1 = pop_int(io,:);
            T(index1) = 1;
            if sum(T) == dim
                Fit(io) = sum(T.*(pr*10).^4)*power(ip);
                Fit_ck(io) = sum(T.*pr)*power(ip);
            else
                Fit(io) = sum(T.*(pr*10).^4)*power(ip) + 100000*abs(dim - sum(T));
                Fit_ck(io) = sum(T.*pr)*power(ip);
            end
            if Fit(io) < best
                best = Fit(io);
                bestz = Fit_ck(io);
                xbest = pop_int(io,:);
                T = zeros(1,96);
                T(xbest) = 1;
                tbest = T;
            end
        end
        record(iter) = bestz;
        FIT1 = 1./Fit;
        sum_Fit = sum(FIT1);
        fitvalue = FIT1./sum_Fit;
        fitvalue = cumsum(fitvalue);
        ms = sort(rand(pop_size, 1));
        fiti = 1;
        newi = 1;
        nf = cell(1, dim);
        while newi <= pop_size
            if(ms(newi) < fitvalue(fiti))
                for im = 1:1:dim
                    nf{1,im}(newi,:) = pop{1,im}(fiti,:);
                end
                newi = newi + 1;
            else
                fiti = fiti + 1;
            end
        end
        for i = 1:2:pop_size
            p = rand;
            if p < Pc
                for im = 1:1:dim
                    q = randi([0,1], 1, chromo_size(im));
                    for j = 1:chromo_size(im)
                        if q(j) == 1
                            temp = nf{1,im}(i+1,j);
                            nf{1,im}(i+1,j) = nf{1,im}(i,j);
                            nf{1,im}(i,j) = temp;
                        end
                    end
                end
            end
        end
        for im = 1:1:dim
            for m = 1:pop_size
                for n = 1:chromo_size(im)
                    r = rand(1,1);
                    if r < Pm
                        nf{1,im}(m,n)=~nf{1,im}(m,n);
                    end
                end
            end
        end
        for im = 1:1:dim
            for m = 1:pop_size
                r1 = rand(1,1);
                if r1 < Pn
                    index = randperm(chromo_size(im),2);
                    temp = nf{1,im}(m,index(1));
                    nf{1,im}(m,index(1)) = nf{1,im}(m,index(2));
                    nf{1,im}(m,index(2)) = temp;
                end
            end
        end
        clear pop
        pop = nf;
    end
elseif kfzd(ip) == 0
    dim = 1;
    upper1 = upper(ip) - t(ip) + 1;
    chromo_size = ceil(log2((upper1 - lower(ip) + 1) / 0.1 + 1));
    pop = cell(1,dim);
    pop_int = zeros(pop_size, dim);
    for iii = 1:dim
        pop{1,iii} = initilize_pop(pop_size, chromo_size);
    end
    record = zeros(1, max_iter);
    best = inf;
    for iter = 1:1:max_iter
        for il = 1:1:dim
            pop_int(:,il) = ...
                round(pop_decode(pop{1,il}, pop_size, chromo_size(il), upper1, lower(ip)));
        end
        Fit = zeros(pop_size, 1);
        Fit_ck = zeros(pop_size, 1);
        for io = 1:1:pop_size
            index1 = pop_int(io,1);
            T = zeros(1,96);
            T(index1 : index1+t(ip)-1) = 1;
            Fit(io) = sum(T.*(10*pr).^4)*power(ip);
            Fit_ck(io) = sum(T.*pr)*power(ip);
            if Fit(io) < best
                best = Fit(io);
                bestz = Fit_ck(io);
                xbest = pop_int(io,:);
                T = zeros(1,96);
                T(xbest : xbest+t(ip)-1) = 1;
                tbest = T;
            end
        end
        record(iter) = bestz;
        FIT1 = 1./Fit;
        sum_Fit = sum(FIT1);
        fitvalue = FIT1./sum_Fit;
        fitvalue = cumsum(fitvalue);
        ms = sort(rand(pop_size, 1));
        fiti = 1;
        newi = 1;
        nf = cell(1, dim);
        while newi <= pop_size
            if(ms(newi) < fitvalue(fiti))
                for im = 1:1:dim
                    nf{1,im}(newi,:) = pop{1,im}(fiti,:);
                end
                newi = newi + 1;
            else
                fiti = fiti + 1;
            end
        end
        for i=1:2:pop_size
            p = rand;
            if p < Pc
                for im = 1:1:dim
                    q = randi([0,1], 1, chromo_size(im));
                    for j = 1:chromo_size(im)
                        if q(j) == 1
                            temp = nf{1,im}(i+1,j);
                            nf{1,im}(i+1,j) = nf{1,im}(i,j);
                            nf{1,im}(i,j) = temp;
                        end
                    end
                end
            end
        end
        for im = 1:1:dim
            for m = 1:pop_size
                for n = 1:chromo_size(im)
                    r = rand(1,1);
                    if r < Pm
                        nf{1,im}(m,n)=~nf{1,im}(m,n);
                    end
                end
            end
        end
        for im = 1:1:dim
            for m = 1:pop_size
                r1 = rand(1,1);
                if r1 < Pn
                    index = randperm(chromo_size(im), 2);
                    temp = nf{1,im}(m, index(1));
                    nf{1,im}(m,index(1)) = nf{1,im}(m,index(2));
                    nf{1,im}(m,index(2)) = temp;
                end
            end
        end
        clear pop
        pop = nf;
    end
end