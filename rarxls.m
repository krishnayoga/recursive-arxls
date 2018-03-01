function [output] = rarxls(data, orde, ff)

	% Menentukan initial value dari covariance matrix P
	alpha = 100;

	% Mendapat ukuran data
	[sizeRow, sizeColumn] = size(data);
	if sizeColumn ~= 2, error('y dan u masing-masing harus terdiri dari satu kolom'); end
	N = sizeRow;

	% Mendefinisikan m sebagai banyak parameter
	% dan I sebagai matriks identitas mxm
	m = orde * 2;
	I = eye(m);

	% Memisahkan y dan u ke vektor masing-masing
	y = data(:, 1);
	u = data(:, 2);

	% Membuat matrix yang nanti digunakan
	phi = zeros(m, 1, N);
	e = zeros(N, 1);
	P = zeros(m, m, N);
	gammaK = zeros(m, 1, N);
	theta = zeros(m, 1, N);
	yhat = zeros(N, 1);
	aperror = zeros(N, 1);
	a = zeros(N, 1);
	b = zeros(N, 1);



	% Melakukan pengulangan dari k=1 sampai k=N
	for k = 1:N

		% Mengisi data vector phi
		for iterOrder = 1:orde
			if k - iterOrder <= 0
				phi(iterOrder, 1, k) = 0;
				phi(iterOrder + orde, 1, k) = 0;
			else
				phi(iterOrder, 1, k) = -1*y(k - iterOrder);
				phi(iterOrder + orde, 1, k) = u(k - iterOrder);
			end
		end

		% Mengisi prediction error e(k)
		if k - 1 <= 0, prevTheta = zeros(m, 1);
		else, prevTheta = theta(:, :, k-1);
		end
		e(k) = y(k) - (phi(:,:,k)' * prevTheta);

		% Mengisi covariance matrix P dan gamma
		if k - 1 <= 0, prevP = alpha * I;
		else, prevP = P(:, :, k-1);
		end
		gammaK(:, :, k) = (prevP * phi(:, :, k)) / (ff + (phi(:, :, k)' * prevP * phi(:, :, k)));
		P(:, :, k) = (1/ff) * (I - (gammaK(:, :, k)*phi(:, :, k)')) * prevP;

		% Mengisi parameter vector theta
		theta(:,:, k) = prevTheta + gammaK(:, :, k) * e(k);

		% Menghitung yhat(k) dan aperror(k)
		yhat(k) = phi(:, :, k)' * theta(:, :, k);
		aperror(k) = y(k) - yhat(k);

		% Mengisi matrix parameter a dan b
		a(k, 1) = theta(1, 1, k);
		b(k, 1) = theta(2, 1, k);
	end



	% Mendefinisikan 'parameter' sebagai struct yang berisi matriks parameter a dan b
	parameter = struct('a', a, 'b', b);

	% Mendefinisikan output dari program sebagai struct
	output = struct('parameter', parameter, 'yhat', yhat, 'aperror', aperror);
end