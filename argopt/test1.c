#include "argopt.h"

int main(int argc, char **argv)
{
  argopt_t *ao;
  const char *fn = NULL;
  char *str = NULL;
  int n = 0, freq = -1, verbose = 0, sw = 0, day = 0;
  double x = 1e-31;
  int nset, xset, fset;
  const char *weekdays[] = {"Sunday", "Monday", "Tuesday",
    "Wednesday", "Thursday", "Friday", "Saturday"};

  ao = argopt_open(0);

  argopt_add(ao, NULL, "!", &fn, "input file");
  /* argopt_add() == argopt_regopt() */
  argopt_add(ao, "--off", "%b", &sw, "turn off");
  argopt_add(ao, "--freq", "%d", &freq, "frequency of saving files");
  argopt_add(ao, "-n", "!%d", &n, "an integer n");
  argopt_add(ao, "-x", "%lf", &x, "a double x");
  argopt_add(ao, "-s", NULL, &str, "a string (allowing , as the separator)");
  argopt_add(ao, "-v", "%+", &verbose, "verbose level, -vv means more verbose");
  argopt_addx(ao, "-d", "%list", &day, "weekday", weekdays, 7);
  argopt_addversion(ao, "--version");
  argopt_addhelp(ao, "-h");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  nset = argopt_isset(ao, n);
  xset = argopt_isset(ao, x);
  fset = argopt_isset(ao, freq);
  printf("fn %s, verbose %d, sw %d, weekday %s, "
      "n[%d] %d, x[%d] %g, freq[%d] %d, str ``%s''\n",
      fn, verbose, sw, weekdays[day], nset,
      n, xset, x, fset, freq, str);
  argopt_close(ao);
  return 0;
}
