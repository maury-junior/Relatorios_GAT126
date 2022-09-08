%% Script para Solução do Trabalho de Identificação
% 2018-1
% Bruno Henrique Groenner Barbosa
% 10/05/2018

%% Limpar Workspace

clear all
close all
clc
warning off

%% Gerar Arquivos de Entrada

% Gerar entrada do tipo degrau
n_amostras = 3000;
degrau = ones(1,n_amostras);

fileID = fopen('degrau.txt','w');
fprintf(fileID,'%6.2f \n',degrau);
fclose(fileID);

% Gerar entrada PRBS
prbs = idinput(n_amostras,'prbs',[0 1],[0 1]);

fileID = fopen('prbs.txt','w');
fprintf(fileID,'%6.2f \n',prbs);
fclose(fileID);

%% Ler dados do Labview
% obs: trocar virgula por ponto antes

dados_degrau = load('ensaio_degrau.txt');
dados_prbs = load('ensaio_prbs.txt');
