clear all 
close all
clc

% EQUACIONAMENTO GERAL DO SISTEMA

% PARÂMETROS DE ENTRADA FORNECIDOS PELO ARTIGO:


% Número de dentes (N)
N_p = 17; % Número de dentes do pinhão
N_c = 26; % Número de dentes da coroa

% Módulo métrico normal (mn)
mn = 4.23;

% raio de filete (r)
raio_filete = 0.2 * mn;

% Ângulo de pressão (theta)
theta = 20;

% Largura da face (b) - Recomendado que esteja entre o intervalo 6*m <= b
% <= 16*m
b_p  = linspace(6 * mn, 16 * mn, 11); % Largura da face do pinhão
% Gera 100 valores igualmente espaçados no intervalo

b_c = b_p; % Largura da face da coroa

% Relação de engrenamento (i)
i = N_c / N_p;

% Torque (T) - 
T_p = 250; % Torque do pinhão
T_c = T_p * i; % Torque da coroa
 
% Razão de Poisson do pinhão e da coroa 
v1 = 0.3;
v2 = 0.3;

% Módulo de elasticidade do pinhão e da coroa 
E1 = 210000;
E2 = 210000;

% Número de contatos por revolução (q)    
q = 1;

% Rotação em rpm (omega) (VERIFICAR)
omega_p  = 8000; % Rotação do pinhão
omega_c = omega_p / i; % Rotação da coroa



% CARACTERÍSTICAS DAS ENGRENAGENS


% Diametro primitivo (dp)
dp_p = 71.91; % Diâmetro primitivo do pinhão
dp_c = 109.98; % Diâmetro primitivo da coroa

% Diâmetro externo (De)
De_p = 80.; % Diâmetro externo da coroa
De_c = 118.44; % Diâmetro externo da coroa

% Diâmetro raiz (Dr)
Dr_p = 61.335; % Diâmetro raiz do pinhão
Dr_c = 99.405; % Diâmetro raiz da coroa

% Diâmetro base (Db)
Db_p = 67.573; % Diâmetro base do pinhão
Db_c = 103.347; % Diâmetro base da coroa

% Distância entre centros (Cdist)
Cdist = (dp_p + dp_c) / 2;

% Passo circular (p)
p = mn * pi;

% Espessura do dente (e)
e = p / 2;

% Velocidade linear (V)
V_p = (pi * dp_p * omega_p)/60000; % Velocidade linear do pinhão
V_c = (pi * dp_c * omega_c)/60000; % Velocidade linear do pinhão

% CALCULAR A POTÊNCIA TRANSMITIDA [W]

% Raio primitivo (rp)
rp_p = dp_p / 2; % Raio primitivo do pinhão
rp_c = dp_c / 2; % Raio primitivo da coroa



% FATORES PARA CÁLCULO SEGUNDO A NORMA AGMA


% Carga tangencial (Ft) - Avaliar qual equação utilizar de acordo com os
% dados de entrada que possui (para alterar eq. verificar a monografia)
Ft = (1000 * T_p) / rp_p; % Carga tangencial transmitida


% Fator de sobrecarga (Ko) - Avaliar o fator de acordo com o tipo de fonte
% de alimentação e máquina de acordo com a Tabela 3 do trabalho.
Ko = 1;

% Fator dinâmico (Kv) - Encontrar via gráfico apresentado na monografia ou
% pela equação abaixo:
Av = 6;
B = 0.25 * ((Av - 5)^(2/3));
C = 3.5637 + 3.9914 *(1 - B);
Kv = (C ./ (C + sqrt(V_p))).^(-B);

% Fator de tamanho (Ks) - Buscar o valor de Y na tabela presente na tese;

% Fator de forma de Lewis (Y)
Y = 0.303;

% Inicializa o vetor Ks com o mesmo tamanho de b_p
Ks = zeros(size(b_p));

% Calcule Ks com base no intervalo de b_p
for i = 1:length(b_p)
    if b_p(i) >= 2.79 && b_p(i) <= 51.99
        Ks(i) = 1 / (1.1833 * (b_p(i) * mn * sqrt(Y))^(-0.0535));
        
    elseif b_p(i) >= 52 && b_p(i) <= 254
        Ks(i) = 1 / (1.4098 * (b_p(i) * mn * sqrt(Y))^(-0.0785));
        
    else
        error('O valor de b_p(%d) está fora do intervalo válido.', i);
    end
end

% Fator de distribuição de carga (Kh)
    % Fator de proporção do pinhão (Khpf)
      
% Inicializa a matriz Khpf com o mesmo tamanho de b_p
Khpf = zeros(size(b_p));

% Percorre cada elemento de b_p
for i = 1:numel(b_p)
    if b_p(i) <= 25
        Khpf(i) = (b_p(i) / (10 * dp_p)) - 0.025;

    elseif b_p(i) > 25 && b_p(i) <= 432
        Khpf(i) = (b_p(i) / (10 * dp_p)) - 0.0375 + 0.000492 * b_p(i);

    else
        error('O valor de b_p(%d) está fora do intervalo válido.', i);
    end
end

     
    % Fator de formato da face do dente (Khmc)
    Khmc = 1; % Sem cooramento no dente
    % Modificador de proporção do pinhão (Khpm)
    Khpm = 1;
    % Fator de alinhamento de engrenagem (Khma)
    A1 = 0.675 * (10^(-1));
    B1 = 0.504 * (10^(-3));
    C1 = -1.44 * (10^(-7));
    Khma = A1 + (B1.*b_p) + (C1 .* (b_p).^2); % Considerando a curva 3 (eng. fechadas de precisão)
    % Fator de ajuste (Khe)
    Khe = 0.8;

Kh = 1.0 + (Khmc * ((Khpf*Khpm) + (Khma*Khe)));

% Fator de confiabilidade (Yz) - Considerando uma confiabilidade de 99%
Yz = 1;   

% Fator de temperatura (Ytheta) - Observar contexto (geralmente mantém 1)
Ytheta = 1;

% Coeficiente de elasticidade (Ze)
    
Ze = sqrt(1 / (pi*((((1 - v1^2) / E1))+((1 - v2^2) / E2))));

% Fator de condição de superfície (Zr) - Considera-se 1
Zr = 1;

% Fator geométrico para resistência ao crateramento (Zi) - Pode ser
% calculado a partir do gráfico mostrado no trabalho ou a partir da equação
% abaixo.
Cc = (i*cosd(theta)*sind(theta)) / (2 * (i +1));
C1 = (N_p * sind(theta) / 2);
C2 = C1 * i;
C3 = pi * cosd(theta); 
C4 = 0.5 * ((sqrt((N_p + 2)^2-(N_p * cosd(theta))^2))-(sqrt((N_p)^2-(N_p * cosd(theta))^2)));
Cx = ((C1 - C3 + C4)*(C2 + C3 - C4)) / (C1 *C2);

Zi = Cc*Cx;

% Fator do ciclo de tensão para resistência ao crateramento (Zn)
    % Buscar valores de L e q
L = 15000;
nl = 60 * L * omega_p * q;
Zn = 1.4488 * (nl^(-0.023));

% Fator de proporção de dureza para resistência ao contato (Zw)
    % Avaliar valor de A na relação apresentada no trabalho
A = 0;    
Zw = 1 + A * (i - 1);

% Fator de espessura de borda (Kb)
    % Espessura de borda abaixo do dente.
tr = 13.14;    
ht = mn + 1.25 * mn;

mb = tr / ht;

if (mb >= 1.2)
    Kb = 1;
else 
    Kb = 1.6 * log(2.242/mb);
end 

% Fator geométrico para resistência à flexão (Yj) - Avaliar gráfico
% presente no trabalho escrito.
Yj = 0.29;

% Fator do ciclo de tensão para resistência à flexão (Yn)
Yn = 1.3558*(nl^(-0.0178));




                    % TENSÕES AO CONTATO E À FLEXÃO


% Tensão ao contato AGMA
sigma_H = Ze.*sqrt((Ft.*Ko.*Kv.*Ks.*Kh.*Zr)./(dp_p.*b_p.*Zi));
% Certifique-se de que todas as variáveis têm dimensões compatíveis


% Tensão à flexão AGMA
sigma_F = (Ft.*Ko.*Kv.*Ks.*Kh.*Kb)./(b_p.*mn.*Yj);





              % NÚMERO DE TENSÃO AO CONTATO E À FLEXÃO PERMITIDO


% Número de tensão ao contato permitido (sigma_HP) - Pode ser obtido pela
% tabela presente no trabalho ou pela equação descrita abaixo.

% Inicialmente considerando um projeto com fator de segurança = 1.

Sh_prov = 1;
Sigma_HP_prov = (sigma_H * Sh_prov * Ytheta * Yz) / (Zn * Zw);


% Calcular as durezas necessárias (estimadas inicialmente)
HB1_contato = (Sigma_HP_prov - 200) / 2.22;
HB2_contato = (Sigma_HP_prov - 237) / 2.41;
HB = 634; % Definido pelo artigo

% Número de tensão ao contato permitido já com a dureza definida -
% Utilizando a Classe 2;

Sigma_HP = (2.41 * HB) + 237;


% Número de tensão à flexão permitido (sigma_FP) - Pode ser obtido pela
% tabela presente no trabalho ou pela equação descrita abaixo.


% Calcular as durezas necessárias (estimadas inicialmente)

Sigma_FP = (0.703 * HB) + 113;



                    % FATOR DE SEGURANÇA AO CONTATO E À FLEXÃO



% Fator de segurança ao contato
Sh = (Zn.*Zw.*Sigma_HP)./(Ytheta.*Yz.*sigma_H);

% Fator de segurança à flexão
Sf = (Yn.*Sigma_FP)./(Ytheta.*Yz.*sigma_F);





% Exibir resultados

format long;  % Ajusta o formato de exibição para longo

% Resultados finais 

fprintf('Valor da tensão ao contato AGMA: %.3f\n', sigma_H);
fprintf('Fator de segurança ao contato AGMA: %.3f\n', Sh);
fprintf('Valor da tensão à flexão AGMA: %.3f\n', sigma_F);
fprintf('Fator de segurança à flexão AGMA: %.3f\n', Sf);

% Plotagem de gráfico

figure(1)
plot(b_p, sigma_H, 'b.')
xlabel('Largura da face (b) [mm] - no intervalo 6*m <= b <= 16*m')
ylabel('Tensão ao contato AGMA σH [MPa]')
title('Tensão ao Contato em Função da Largura da Face')

% Encontrar o índice do valor próximo a 50.76 em b_p
[~, idx] = min(abs(b_p - 50.76));

% Destacar o ponto onde b_p é aproximadamente 50.76
hold on
plot(b_p(idx), sigma_H(idx), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5)

% Adicionar uma anotação para o ponto específico
text(b_p(idx), sigma_H(idx), sprintf('b_p = %.2f, \\sigma_H = %.2f MPa', b_p(idx), sigma_H(idx)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'red');
hold off




figure(2)
plot(b_p, Sh, 'b.')
xlabel('Largura da face (b) [mm] - no intervalo 6*m <= b <= 16*m')
ylabel('Fator de Segurança ao Contato (Sh)')
title('Fator de Segurança ao Contato em Função da Largura da Face')

% Encontrar o índice do valor próximo a 50.76 em b_p
[~, idx] = min(abs(b_p - 50.76));

% Destacar o ponto onde b_p é aproximadamente 50.76
hold on
plot(b_p(idx), Sh(idx), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5)

% Adicionar uma anotação para o ponto específico
text(b_p(idx), Sh(idx), sprintf('b_p = %.2f, Sh = %.2f', b_p(idx), Sh(idx)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'red');
hold off




figure(3)
plot(b_p, sigma_F, 'b.')
xlabel('Largura da face (b) [mm] - no intervalo 6*m <= b <= 16*m')
ylabel('Tensão de Flexão AGMA σF [MPa]')
title('Tensão à Flexão em Função da Largura da Face')

% Encontrar o índice do valor próximo a 50.76 em b_p
[~, idx] = min(abs(b_p - 50.76));

% Destacar o ponto onde b_p é aproximadamente 50.76
hold on
plot(b_p(idx), sigma_F(idx), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5)

% Adicionar uma anotação para o ponto específico
text(b_p(idx), sigma_F(idx), sprintf('b_p = %.2f, \\sigma_F = %.2f MPa', b_p(idx), sigma_F(idx)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'red');
hold off



figure(4)
plot(b_p, Sf, 'b.')
xlabel('Largura da face (b) [mm] - no intervalo 6*m <= b <= 16*m')
ylabel('Fator de Segurança à Flexão (Sf)')
title('Fator de Segurança à Flexão em Função da Largura da Face')

% Encontrar o índice do valor próximo a 50.76 em b_p
[~, idx] = min(abs(b_p - 50.76));

% Destacar o ponto onde b_p é aproximadamente 50.76
hold on
plot(b_p(idx), Sf(idx), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5)

% Adicionar uma anotação para o ponto específico
text(b_p(idx), Sf(idx), sprintf('b_p = %.2f, Sf = %.2f', b_p(idx), Sf(idx)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'red');
hold off
