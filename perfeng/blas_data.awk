BEGIN{first=0;ind=0;}
{
  data[ind] = ($2 + $5 + $8)/3;
  ind++;
}
END{
  sum=0;
  printf("rollout = [\n");
  for (i=0; i<8; i++) {
    printf("%s\n", data[i]);
    sum+=data[i];
  }
  print "];"
  printf("tot-rollout = %s;\n\n", sum);
}
