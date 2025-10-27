FindRoots::ResultType FindRoots::find0(FunctionType f,double a, double b,double *x)
{
  double Fp,Fm,Fh, step=Step;

  //нахождение промежутка в котором есть корень

  //пропуск промежутка где функция не имет значений
  *x=a-Step;
  do{
    *x+=Step;
    Fp=(FuncClass->*f)(*x);
    if(*x>b) return NOTEXIST;
  }while(!_finite(Fp));
  
  //поиск точной границы промежутка несуществования функции
  if(*x>a+1.0e-14)
  {
    while(step>XAccur){
      step*=0.5;
      *x-=step;
      if(!_finite((FuncClass->*f)(*x)))
        *x+=step;
    }
    step=Step;
    Fp=(FuncClass->*f)(*x);
  }

  if(fabs(Fp)<1.0e-14)
  {
    return EXIST;
  }

  //нахождение промежутка знакопеременности функции
  do{
    Fm=Fp;
    *x+=step;
    if(*x>b)
    {
      step=b-(*x-step);
      *x=b;
    }
    Fp=(FuncClass->*f)(*x);
  }while((Fm*Fp>0.0)&&((*x)<b)&&(_finite(Fp)));

  //поиск точной границы промежутка где функция не имет значений
  if(!_finite(Fp)) 
  {
    *x-=step;
    double tmp=*x;
    do{
      step*=0.5;
      *x+=step;
      if(!_finite((FuncClass->*f)(*x)))
        *x-=step;
    }while(step>XAccur);
    Fp=(FuncClass->*f)(*x);

    if(fabs(Fp)<1.0e-14)
      return EXIST;
    if(Fp*Fm>0.0)
    {
      *x+=Step;  //вообщето должно быть step, но это не работает при точностях <10e-14
      return NOTEXIST;
    }
    step=*x-tmp;
  }

  
  //если промежуток найден, то далее ищется корень (делением пополам)
  if(Fm*Fp<=0.)
  {
//    *x-=step;
    double Fm_Fp;
    while((step>XAccur)&&(Fp!=0.0))
    {
      Fm_Fp=fabs(Fm-Fp);
      step*=0.5;
      *x=*x-step;
      Fh=(FuncClass->*f)(*x);
      if(!_finite(Fh)) return MATHERROR;
      if(Fh*Fm>0.)
      {
        Fm=Fh;
        *x=*x+step;
      }
      else
        Fp=Fh;
    }
    if((fabs(Fp-Fm)<YAccur) || (Fp==0.0))
      return EXIST; //корень есть
    else
    {
      if(Fm_Fp>fabs(Fm-Fp))
        return EXIST;
      return SINGULAR; //разрыв
    }
  }
  else
    return NOTEXIST; // корней нет
}