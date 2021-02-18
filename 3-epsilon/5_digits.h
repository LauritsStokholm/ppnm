 /* EXERCISE 4: OPTIONAL */
void name_digit(int i){

  switch (i)
  {
    case 0:
      printf("zero\n");
      break;
    case 1:
      printf("one\n");
      break;
    case 2:
      printf("two\n");
      break;
    case 3:
      printf("three\n");
      break;
    case 4:
      printf("four\n");
      break;
    case 5:
      printf("five\n");
      break;
    case 6:
      printf("six\n");
      break;
    case 7:
      printf("seven\n");
      break;
    case 8:
      printf("eigth\n");
      break;
    case 9:
      printf("nine\n");
      break;
    case 10:
      printf("ten\n");
      break;

    default:
      printf("not a digit\n");
  }
}

int digits(void){
  for(int i=0; i<11; i++){name_digit(i);}
  return 0;
}
