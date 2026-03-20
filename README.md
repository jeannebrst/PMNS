# Toeplitz with NEON (sizes 14, 15, 16)

**--------------------- WORK STILL IN PROGRESS ---------------------**

This project implements **modular multiplication** using **Toeplitz matrices** optimized with **ARM NEON instructions** for sizes:

- `n = 14`
- `n = 15`
- `n = 16`

This project has been tested on a **Raspberry Pi 400**.

## 1. Requirements

### 1.1 MPHELL library

You need the MPHELL library:

https://www-fourier.univ-grenoble-alpes.fr/~mphell/

### 1.2 AMNS Generator

You also need the AMNS generator (to be installed on a machine where SageMath is available):

https://github.com/arithPMNS/generalisation_amns

Edit the function `get_amns_params()` in `util.c` so that it **calls the generator** and returns parameters in the **required output format** below :

Example :
```c
"[15, 16, [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], 18, 1732996756754789471339355312096369946981361386195029723648, [-2239, -494, -130, -54, -772, 615, -474, -1531, 1849, -895, -1125, -289, 438, -1978, 590, 1077], [4092435061, 3150241284, 3437500926, 2483521910, 1913883334, 2800600421, 2868738738, 3215995195, 1101386429, 229037773, 2613946518, 4042815467, 3393698058, 1925116940, 1693409253, 1245900174], [-5701, -3453, -5708, -4336, -4736, -3772, -3384, -5334, -4480, -2794, -2698, -4467, -761, -3023, -3378, -1571], [-1812, -10495, 232, -4276, -8242, -1249, -1563, -6222, -3852, 97, -6506, -2222, 120, -4246, -3666, -1254]]"
```

It is mandatory that the output of `get_amns_param()` matches the expected structure, even if you need to modify the generator itself.

A reference file to get the generator to match the expected output is given in this repository:

```
amns_generator.py
```

This file is **not meant to be executed**.  
It only serves as a **template/reference for the expected output format**.

## 2. Usage

### 2.1 Modular multiplication

To compute:

```
a × b mod p
```

for a given size `{size} ∈ {14, 15, 16}`:

```bash
make {size}
```

Then enter your values:

```
a b p ? (or 'd' for default values)
```

- Type values manually: `a b p`
- Or type `d` to use default values (editable in `mult_{size}/main_{size}.c`)

### 3.2 Example

Compute:

```
10 × 10 mod 115792089210356248762697446949407573530086143415290314195533631308867097853951
```

with `n = 14`:

```bash
make 14
```

Input:

```
10 10 115792089210356248762697446949407573530086143415290314195533631308867097853951
```

Output:

```
a = 10
b = 10

Normal : a*b = 100
NEON   : a*b = 100
```

## 4. Performance testing

To benchmark performance with random `a` and `b`:

```bash
make perfs_{size}
```

- `{size}` ∈ `{14, 15, 16}`
- `p` is defined in:
  ```
  mult_{size}/main_{size}_perfs.c
  ```

## 5. Notes

- NEON implementation is designed to **accelerate Toeplitz-based multiplication**
- Schoolbook multiplication is still used as a reference
- Results from both implementations are compared for validation
- Users must generate AMNS parameters externally and integrate them manually
