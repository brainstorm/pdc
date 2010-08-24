BEGIN{first=0;ind=0;}
/BLOCK/{
  if (first>0) {
    printf("block%s-%s = [\n", block1, block2);
    sum=0;
    for (i=0; i<8; i++) {
      printf("%s\n", data[i]);
      sum+=data[i];
    }
    print "];"
    printf("tot%s-%s = %s;\n\n", block1, block2, sum);
    ind=0;
  } else {
    first = 1;
  }
  block1=$3;
  block2=$6;
}
!/BLOCK/{
  data[ind] = ($2 + $5 + $8)/3;
  ind++;
}
END{
  printf("block%s-%s = [\n", block1, block2);
  sum=0;
  for (i=0; i<8; i++) {
    printf("%s\n", data[i]);
    sum+=data[i];
  }
  print "];"
  printf("tot%s-%s = %s;\n\n", block1, block2, sum);
}
