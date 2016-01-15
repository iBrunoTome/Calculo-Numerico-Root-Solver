/*
    Nome: root-solver.bin
    Copyright (c) 2014 IFMG. All rights reserved.

    =================== Membros ===================

    Bruno Tomé - 0011254 - ibrunotome@gmail.com
    Matheus Calixto - 0011233 - calixtinn@gmail.com

    =========== Instruções de compilação ==========

    Abra o Terminal e digite:

    cd <DIRETÓRIO>
    gcc root-solver.c -oroot-solver.bin
    ./root-solver.bin <ARQUIVO ENTRADA.TXT> <ARQUIVO SAIDA.HTML>

    ========== Ambiente de Desenvolvimento ========

    Bruno Tomé

    Sistema Operacional: OS X 10.9.4
    Hardware: Core i5 3ª Geração 2.5 GHz| 16 GB RAM 1600Mhz
    Desenvolvido no Xcode 5.1.1 e no editor de textos Sublime Text 3
    Compilado no Terminal do OS X

    GCC

    Apple LLVM version 5.1 (clang-503.0.40) (based on LLVM 3.4svn)
    Target: x86_64-apple-darwin13.3.0
    Thread model: posix

    Matheus Calixto

    Sistema Operacional: Linux Ubuntu 14.04 LTS
    Hardware: Core i5 3ª Geração 1.8 ~ 2.5 GHz| 8GB RAM 1600 Mhz
    Desenvolvido no NetBeans 8.0
    Compilado no Terminal do Ubuntu.

    GCC (Ubuntu 4.8.2-19ubuntu1) 4.8.2

    ============ Objetivo do programa ============

    Resolver equações polinomiais pelos métodos de Newton, Bissecção e Regula Falsi
    e exibir os resultados em um arquivo html.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int reserva = 0, flag = 0, expoente = 0, iterations = 0, poligrau = 0, cont = 0, temp;
float limdown = 0, limtop = 0, error = 0, root = 0;
char expression[256], coeftemp[256];
float vetorfunc[256], vetorderiv[256],raiz;

float potencia(float a, int b){
    float pot = 1;
    int i = 0;
    
    for (i = 0; i < b; i++){
        pot *= a;
    }
    
    return pot;
}

float polinomio(float x){
    int i = 0;
    float soma = 0;
    expoente = reserva;
    
    while (expoente != -1){
        soma = (vetorfunc[i] * potencia(x,expoente) + soma);
        i++;
        expoente--; 
    }
    return soma;
}

// Funções metodohorner, bisseccao, regulaFalsi e newthonraphson implementadas em C a partir do pseudo-código dos slides

// Implementação Método Horner
float metodohorner(float c[],float a){
    float y;
    int i = 1, n;
    n = poligrau;
    y = c[0];
    
    do{
        y =  y * a + c[i];
        i++;
    }while(i<(n+1));
    
    return y; 
}

float metodohornerderiv(float c[],float a){
    float y;
    int i = 1, n;
    n = poligrau-1;
    y = c[0];
    
    do{
        y =  y * a + c[i];
        i++;
    
    }while(i<(n + 1));
    
    return y;
}

// Fim Método Horner

// Implementação Método da Bissecção

float bisseccao(float a, float b, float epsilon, int maxiter, FILE *html){
    float delta = 0, x = 0,xa,xb,xx;
    int condicao, iter = 0;   
    
    if ((polinomio(a) * polinomio(b)) > 0){
        fprintf(html,"<h3>Função não muda o sinal - Erro</h3>");
        flag=1;
    }
    else flag=0;
    
    delta = (b - a) / 2;
    iter = 0;
    condicao = 1;
    
    while(condicao == 1) {
        
        x = (a + b) / 2;
        xa=polinomio(a);
        xb=polinomio(b);
        xx=polinomio(x);
        
        fprintf(html,"<tr align=center><td>%d</td><td>%f</td><td>%f</td><td>%f</td><td>%f</td><td>%f</td><td>%f</td><td>%f</td></tr>",iter,a,xa,b,xb,x,xx,delta);
        
        if ((delta <= epsilon)  && (polinomio(x) < epsilon)){
            condicao = 0;
        }
        
        else if (iter > maxiter) {
            condicao=0;
        }
        
        else if ((polinomio(a) * polinomio(x)) > 0) {
            a = x;
        }
        else {
            b = x;
        }
        delta = delta / 2;
        iter++;
    }
    return x; 
}

// Fim Método Bissecção

// Implementação Método Regula Falsi

float regulaFalsi(float b, float a, float epsilon, int maxiter,FILE *html){
    
    float x = 0, aux = 0, delta = 0;
    int iter,condicao;
    
    if ((polinomio(a) * polinomio(b)) > 0){
        fprintf(html,"<h3>Função não muda o sinal - Erro</h3>");
        flag=1;
    }
    
    else  flag = 0; 
    
    if (polinomio(a) > 0){
        aux = a;
        a = b;
        b = aux;
    }
    
    iter = 0;
    x = b;
    condicao=1;
    
    while (condicao != 0){
        
        delta = -((polinomio(x) / (polinomio(b) - polinomio(a)) * (b - a)));
        
        x = x + delta;
        
        fprintf(html,"<tr align=center><td>%d</td><td>%f</td><td>%f</td><td>%f</td><td>%f</td><td>%f</td><td>%f</td><td>%f</td></tr>",iter,a,polinomio(a),b,polinomio(b),x,polinomio(x),delta);
        
        if (iter >= maxiter){
            condicao = 0;
        }
        
        if (polinomio(x) < 0){
            a = x;
        }
        else {
            b = x;
        }
        
        iter++;
        
        if  ((fabs(delta) <= epsilon) && (fabs(polinomio(x)) <= epsilon)) {
            condicao=0;
        }       
    }
    
    return x;
}

// Fim Método Regula Falsi

// Implementação Método Newton Raphson

float newtonraphson(float x0, float epsilon, int maxiter,FILE *html, float c[],float d[]){
    int iter = 0, repetir = 0;
    float x = x0, fx, dx, delta = 0;
    
    fx=metodohorner(c,x0);
    dx=metodohornerderiv(d,x0);
    
    fprintf(html,"<tr align=center><td>%d</td><td>%f</td><td>%f</td><td>%f</td><td> ----- </td></tr>",iter,x,fx,dx);
    repetir = 1;
    flag = 0;
    while (repetir == 1){
        delta = -fx/dx;
        x += delta;
        fx = metodohorner(c,x);
        dx = metodohornerderiv(d,x);
        iter++;
        
        fprintf(html,"<tr align=center><td>%d</td><td>%f</td><td>%f</td><td>%f</td><td>%f</td></tr>",iter,x,fx,dx,delta);
        
        if ((dx == 0) || (iter >= maxiter )){
            repetir = 0;
            fprintf(html,"Método de Newton nao convergiu\n");
            fprintf(html,"<p><h3> Método de Newton nao convergiu, valor de x não é preciso </h3></p>");
            flag = 1;
        }
        
        if ((fabs(delta) <= epsilon) && (fabs(fx) <= epsilon)){
            repetir = 0;
        }
    }
    
    return x;
}

// Fim Método Newton Raphson

// Lerndo arquivo de entrada

void lendoarquivo(char **argv) {
    
    int i=0, k = 0, tamanho = 0,auxtam=0,j=0,u,tamanho2=0;
    char *veriflinha, auxcoef[256], linha[256],iteracoes[256], inferior[256],superior[256],erro[256],grau[256],valorcoef[256];
    
    FILE *arquivoin;
    
    arquivoin = fopen(argv[1],"r");
    
    while (!feof(arquivoin)){
        veriflinha = fgets (linha, 256, arquivoin);
        if (veriflinha) {
            
            if (linha[0] == 'd'){
                tamanho=strlen(linha);
                k=0;
                
                for (i=2;i<(tamanho-1);i++) {
                    expression[k]=linha[i];
                    k++;
                    
                } // end for
                
            } // end if linha = d            
            
            else if (linha[0] == 'm') {
                tamanho = strlen(linha);
                k = 0;
                for(i=2;i<(tamanho-1);i++) {
                    iteracoes[k]=linha[i];
                    k++;
                } //end for verificação M
                
            } // end if linha = m
            
            else if (linha[0] == 'l') {
                tamanho = strlen(linha);
                k = 0;
                for(i=2;i<(tamanho-1);i++) {
                    inferior[k]=linha[i];
                    k++;
                } // end for verificação I
                
            } // end if linha = i
            
            else if (linha[0] == 'u') {
                tamanho = strlen(linha);
                k = 0;
                for(i=2;i<(tamanho-1);i++) {
                    superior[k]=linha[i];
                    k++;
                } //end for verificação U
                
            } // end if linha = u
            
            else if (linha[0] == 'e') {
                tamanho = strlen(linha);
                k = 0;
                for(i=2;i<(tamanho-1);i++) {
                    erro[k]=linha[i];
                    k++;
                } //end for verificação E
                
            } // end if linha = e
            
            else if (linha[0] == 'n') {
                tamanho = strlen(linha);
                k = 0;
                for(i=2;i<(tamanho-1);i++) {
                    grau[k]=linha[i];
                    k++;
                } //end for verificação n
                
            } // end if linha = n
                  
            else if (linha[0] == 'a') {
                tamanho = strlen(linha);
                k = 0;
                for(i=2;i<(tamanho-1);i++) {
                    valorcoef[k]=linha[i];
                    k++;
                } //end for verificação A
                
            } // end if linha = a   
            
        } // end verif linha
        
    } // end while
    
    fclose(arquivoin); // Fechando o arquivo de Leitura.
    
    iterations=atoi(iteracoes); // Salva o Número Máximo de Iterações em uma variável inteira...
    poligrau=atoi(grau); // Salva o Grau do Polinômio em uma variável inteira...
    limdown=atof(inferior); // Salva o Limite Inferior em uma variável do tipo Float...
    limtop=atof(superior); // Salva o Limite Superior em uma variável do tipo Float...
    error=atof(erro); // Salva o Erro em uma variável do tipo Float...
    reserva = poligrau; // Salva o grau do polinômio.
    valorcoef[tamanho-3] = valorcoef[tamanho - 3] + ' '; // Adiciona um arroba no final da string para controle.
    tamanho=strlen(valorcoef);
    temp=poligrau;
       
    cont = 0;
    u = 0;
    i = 0;
    
    // Verificações arquivo de entrada

    while (valorcoef[i] != '\0') {
        
        if (valorcoef[i] != ' '){
            coeftemp[u]=valorcoef[i];
            i++;
            u++;
        }
        
        else {
            vetorfunc[cont]=atof(coeftemp);
            i++;
            u = 0;
            cont++;
            strcpy(coeftemp,"");
        }
    }

    for(i = 0;i<poligrau;i++) {
        vetorderiv[i]=temp*vetorfunc[i];
        temp--;
    } 
    
} //end função ler arquivo...

int main(int argc, char **argv){
    
    lendoarquivo(argv);
    
    FILE *html;
    
    html=fopen(argv[2],"w+");
    
    // Impressão HTML
    fprintf(html,"<head><title>Trabalho prático I</title><meta http-equiv=Content-Type content=text/html;charset=utf-8 ><style type=text/css>h2{font-weight: 100}h3{font-weight: 100}table, td, tr{border: solid 1px #000;}</style></head>");
    fprintf(html,"<body>");
    fprintf(html,"<h3>Arquivo de Saída: saida-solver.html</h3>");
    fprintf(html,"<hr align=center size=1>");
    fprintf(html,"<h2 align=center><i>%s</i></h1>",expression);
    fprintf(html,"<hr align=center size=2>");
    
    // Imprime Bisecção
    fprintf(html,"<h2><i>Método: Bissecção</i></h2>");
    fprintf(html,"<table align=center border=1 cellspacing=0 width=80%%>");
    fprintf(html,"<tr align=center><td> <b>Iter</b></td><td><b>l</b></td><td><b>f(l)</b></td><td><b>u</b></td><td><b>f(u)</b></td><td><b>x</b></td><td><b>f(x)</b></td><td><b>Delta</b></td></tr>");
    raiz = bisseccao(limdown,limtop,error,iterations,html);
    fprintf(html,"</table>");
    if (flag == 0)fprintf(html,"<p><h3> Encontrou a raíz real %f com precisão epsilon = %.3f </h3></p>",raiz,error);
    
    // Imprime RegulaFalsi
    fprintf(html,"<h2><i>Método: Regula Falsi</i></h2>");
    fprintf(html,"<table align=center border=1 cellspacing=0 width=80%%/>");
    fprintf(html,"<tr align=center><td><b>Iter</b></td><td><b>l</b></td><td><b>f(l)</b></td><td><b>u</b></td><td><b>f(u)</b></td><td><b>x</b></td><td><b>f(x)</b></td><td><b>Delta</b></td></tr>");
    raiz = regulaFalsi(limtop, limdown, error, iterations,html);
    fprintf(html,"</table>");
    if (flag == 0)fprintf(html,"<p><h3> Encontrou a raíz real %f com precisão epsilon = %.3f </h3></p>",raiz,error);
    
    // Imprime Newton-Raphson
    fprintf(html,"<h2><i>Método: Newton-Raphson </i></h2>");
    fprintf(html,"<table align=center border=1 cellspacing=0 width=50%%>");
    fprintf(html,"<tr align=center><td><b>Iter</b></td><td><b>x</b></td><td><b>f(x)</b></td><td><b>f'(x)</b></td><td><b>Delta</b></td></tr>");
    raiz = newtonraphson(limdown, error, iterations,html,vetorfunc,vetorderiv);
    fprintf(html,"</table>");
    if (flag == 0)fprintf(html,"<p><h3>Encontrou a raíz real %f com precisão epsilon = %.3f </h3></p>",raiz,error);
    fprintf(html,"<h3><center><b><i>Arquivo gerado automaticamente pelo aplicativo</i></b> <i>root-solver.bin</i></center></h3>");
    fprintf(html,"</body>");
    fprintf(html,"</html>");
    
    fclose(html);
    return 0;
}