function [activeX, activeY, activeZ] = sampler_update_flipswap(B, activeX, activeY, activeZ, D, params, outputflag)

    %% emitter sellection
    pairs = get_all_pairs(B);

    displ_sqrd = diff(activeX).^2 + diff(activeY).^2 + diff(activeZ).^2;
    diff_in_t = diff(params.t_mid(:));

    %% samplers
    for i = randperm(size(pairs, 1))

        for n = randperm(params.N)

            for k = randperm(params.K)

                curr_t = params.t_idx(n, k);

                if curr_t == 1
                    continue
                else
                    prev_t = curr_t - 1;
                end

                log_r = displ_sqrd(prev_t, pairs(i, 1)) + displ_sqrd(prev_t, pairs(i, 2));

                [Xprop, Yprop, Zprop, tag] = get_traj_prop( ...
                    activeX(curr_t:end, [pairs(i, 1) pairs(i, 2)]), ...
                    activeY(curr_t:end, [pairs(i, 1) pairs(i, 2)]), ...
                    activeZ(curr_t:end, [pairs(i, 1) pairs(i, 2)]));

                log_r = log_r ...
                    - (Xprop(1, 1) - activeX(prev_t, pairs(i, 1)))^2 ...
                    - (Yprop(1, 1) - activeY(prev_t, pairs(i, 1)))^2 ...
                    - (Zprop(1, 1) - activeZ(prev_t, pairs(i, 1)))^2 ...
                    - (Xprop(1, 2) - activeX(prev_t, pairs(i, 2)))^2 ...
                    - (Yprop(1, 2) - activeY(prev_t, pairs(i, 2)))^2 ...
                    - (Zprop(1, 2) - activeZ(prev_t, pairs(i, 2)))^2;

                log_r = log_r / (4 * D * diff_in_t(prev_t));

                if log(rand) < log_r

                    if outputflag

                        if tag(1)
                            disp([ ...
                                    '---> flip: at n=', num2str(n), ...
                                    ' and k=', num2str(k), ...
                                    ', emitter ', num2str(pairs(i, 1)) ...
                                ])
                        end

                        if tag(2)
                            disp([ ...
                                    '---> flip: at n=', num2str(n), ...
                                    ' and k=', num2str(k), ...
                                    ', emitter ', num2str(pairs(i, 2)) ...
                                ])
                        end

                        if tag(3)
                            disp([ ...
                                    '---> swap: at n=', num2str(n), ...
                                    ' and k=', num2str(k), ...
                                    ', between ', num2str(pairs(i, 1)), ...
                                    ' and ', num2str(pairs(i, 2)) ...
                                ])
                        end

                    end

                    activeX(curr_t:end, [pairs(i, 1) pairs(i, 2)]) = Xprop;
                    activeY(curr_t:end, [pairs(i, 1) pairs(i, 2)]) = Yprop;
                    activeZ(curr_t:end, [pairs(i, 1) pairs(i, 2)]) = Zprop;
                    displ_sqrd = diff(activeX).^2 + diff(activeY).^2 + diff(activeZ).^2;

                end

            end

        end

    end

end

function pairs = get_all_pairs(N)
    pairs = zeros(N * (N - 1) / 2, 2);
    [elem1, elem2] = ndgrid(1:N);

    elem1 = tril(elem1, -1);
    elem2 = tril(elem2, -1);

    pairs(:, 1) = elem1(elem1 > 0);
    pairs(:, 2) = elem2(elem2 > 0);
end

function [X, Y, Z, tag] = get_traj_prop(X, Y, Z)

    tag = false(1, 3);

    switch randi(4)
        case 2
            Z(:, 1) = -Z(:, 1);
            tag(1) = true;
        case 3
            Z(:, 2) = -Z(:, 2);
            tag(2) = true;
        case 4
            Z(:, 1) = -Z(:, 1);
            Z(:, 2) = -Z(:, 2);
            tag(1:2) = true;
    end

    if rand > 0.5
        X(:, [1 2]) = X(:, [2 1]);
        Y(:, [1 2]) = Y(:, [2 1]);
        Z(:, [1 2]) = Z(:, [2 1]);
        tag(3) = true;
    end

end
