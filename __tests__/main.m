L = 1;
dx = L / 10;

dxs = [dx / 2, dx, dx, dx, dx, dx, dx, dx, dx, dx, dx / 2];
mis = [1, 1, 1.2, 1.2, 1.4, 1.4, 1.2, 1.2, 1, 1] / 10;
EIs = [1, 1, 1.22, 1.44, 1.70, 1.96, 1.7, 1.44, 1.22, 1, 1];

[v, wn] = myklestad_pinned_pinned(mis, EIs, dxs);
