//
.att_syntax
.text

.p2align 4
.global secp256k1_fe_mul_inner
.type	secp256k1_fe_mul_inner, %function
secp256k1_fe_mul_inner:
movq	%rsp, %r11
subq	$64, %rsp

movq	%r11,  0(%rsp)
movq	%r12,  8(%rsp)
movq	%r13, 16(%rsp)
movq	%r14, 24(%rsp)
movq	%r15, 32(%rsp)
movq	%rbp, 40(%rsp)
movq	%rbx, 48(%rsp)
movq	%rdi, 56(%rsp)

movq	%rdx, %rcx

movq 	8(%rsi),%rax
mulq 	32(%rcx)
movq 	%rax, %r8
movq 	%rdx, %rbx

movq 	16(%rsi),%rax
mulq 	24(%rcx)
addq 	%rax, %r8
adcq 	%rdx, %rbx

movq 	24(%rsi),%rax
mulq 	16(%rcx)
addq 	%rax, %r8
adcq 	%rdx, %rbx

movq 	32(%rsi),%rax
mulq 	8(%rcx)
addq 	%rax, %r8
adcq 	%rdx, %rbx

movq 	16(%rsi),%rax
mulq 	32(%rcx)
movq 	%rax, %r10
movq 	%rdx, %r11

movq 	24(%rsi),%rax
mulq 	24(%rcx)
addq 	%rax, %r10
adcq 	%rdx, %r11

movq 	32(%rsi),%rax
mulq 	16(%rcx)
addq 	%rax, %r10
adcq 	%rdx, %r11

movq 	24(%rsi),%rax
mulq 	32(%rcx)
movq 	%rax, %r12
movq 	%rdx, %r13

movq 	32(%rsi),%rax
mulq 	24(%rcx)
addq 	%rax, %r12
adcq 	%rdx, %r13

movq 	32(%rsi),%rax
mulq 	32(%rcx)
movq 	%rax, %r14
movq 	%rdx, %r15

shld    $12, %r8, %rbx
shld    $12, %r10, %r11
shld    $12, %r12, %r13
shld    $12, %r14, %r15

movq    $0xFFFFFFFFFFFFF, %rbp

andq	%rbp, %r8
andq	%rbp, %r10
andq	%rbp, %r12
andq	%rbp, %r14

addq    %rbx, %r10
addq    %r11, %r12
addq    %r13, %r14

movq    $0x1000003D10, %rdi

movq    %r8, %rax
mulq    %rdi
movq    %rax, %r8
movq    %rdx, %r9

movq 	0(%rsi),%rax
mulq 	0(%rcx)
addq 	%rax, %r8
adcq 	%rdx, %r9

shld    $12, %r8, %r9

movq    %r10, %rax
mulq    %rdi
movq    %rax, %r10
movq    %rdx, %r11

movq 	0(%rsi),%rax
mulq 	8(%rcx)
addq 	%rax, %r10
adcq 	%rdx, %r11

movq 	8(%rsi),%rax
mulq 	0(%rcx)
addq 	%rax, %r10
adcq 	%rdx, %r11

shld    $12, %r10, %r11

movq    %r12, %rax
mulq    %rdi
movq    %rax, %r12
movq    %rdx, %r13

movq 	0(%rsi),%rax
mulq 	16(%rcx)
addq 	%rax, %r12
adcq 	%rdx, %r13

movq 	8(%rsi),%rax
mulq 	8(%rcx)
addq 	%rax, %r12
adcq 	%rdx, %r13

movq 	16(%rsi),%rax
mulq 	0(%rcx)
addq 	%rax, %r12
adcq 	%rdx, %r13

shld    $12, %r12, %r13

movq    %r14, %rax
mulq    %rdi
movq    %rax, %rbx
movq    %rdx, %rbp

movq 	0(%rsi),%rax
mulq 	24(%rcx)
addq 	%rax, %rbx
adcq 	%rdx, %rbp

movq 	8(%rsi),%rax
mulq 	16(%rcx)
addq 	%rax, %rbx
adcq 	%rdx, %rbp

movq 	16(%rsi),%rax
mulq 	8(%rcx)
addq 	%rax, %rbx
adcq 	%rdx, %rbp

movq 	24(%rsi),%rax
mulq 	0(%rcx)
addq 	%rax, %rbx
adcq 	%rdx, %rbp

shld    $12, %rbx, %rbp

movq    %r15, %rax
mulq    %rdi
movq    %rax, %r14
movq    %rdx, %r15

movq 	0(%rsi),%rax
mulq 	32(%rcx)
addq 	%rax, %r14
adcq 	%rdx, %r15

movq 	8(%rsi),%rax
mulq 	24(%rcx)
addq 	%rax, %r14
adcq 	%rdx, %r15

movq 	16(%rsi),%rax
mulq 	16(%rcx)
addq 	%rax, %r14
adcq 	%rdx, %r15

movq 	24(%rsi),%rax
mulq 	8(%rcx)
addq 	%rax, %r14
adcq 	%rdx, %r15

movq 	32(%rsi),%rax
mulq 	0(%rcx)
addq 	%rax, %r14
adcq 	%rdx, %r15

shld    $16, %r14, %r15

movq    $0xFFFFFFFFFFFFF, %rsi
movq    $0xFFFFFFFFFFFF, %rcx
movq    $0x1000003D1, %rdi

andq	%rsi, %r8
andq	%rsi, %r10
andq	%rsi, %r12
andq	%rsi, %rbx
andq	%rcx, %r14

addq    %r9, %r10
addq    %r11, %r12
addq    %r13, %rbx
addq    %rbp, %r14

movq    %rbx, %rdx
shrq    $52, %rdx
addq    %rdx, %r14
andq	%rsi, %rbx

movq    %r14, %rax
shrq    $48, %rax
addq    %r15, %rax
andq	%rcx, %r14

mulq    %rdi
addq    %rax, %r8
adcq    $0, %rdx
shld    $12, %r8, %rdx
andq	%rsi, %r8
addq    %r10, %rdx

movq    %rdx, %r10
shrq    $52, %rdx
addq    %r12, %rdx
andq	%rsi, %r10

movq    %rdx, %r12
shrq    $52, %rdx
addq    %rbx, %rdx
andq	%rsi, %r12

movq    %rdx, %r13
shrq    $52, %rdx
addq    %rdx, %r14
andq	%rsi, %r13

movq 	56(%rsp), %rdi

movq 	%r8,   0(%rdi)
movq 	%r10,  8(%rdi)
movq 	%r12, 16(%rdi)
movq 	%r13, 24(%rdi)
movq 	%r14, 32(%rdi)

movq	 8(%rsp), %r12
movq	16(%rsp), %r13
movq	24(%rsp), %r14
movq	32(%rsp), %r15
movq	40(%rsp), %rbp
movq	48(%rsp), %rbx
movq	 0(%rsp), %rsp
ret
.size secp256k1_fe_mul_inner, .-secp256k1_fe_mul_inner

.p2align 4
.global secp256k1_fe_sqr_inner
.type secp256k1_fe_sqr_inner, %function
secp256k1_fe_sqr_inner:
movq    %rsp, %r11
subq    $80, %rsp

movq    %r11,  0(%rsp)
movq    %r12,  8(%rsp)
movq    %r13, 16(%rsp)
movq    %r14, 24(%rsp)
movq    %r15, 32(%rsp)
movq    %rbp, 40(%rsp)
movq    %rbx, 48(%rsp)
movq    %rdi, 56(%rsp)

movq 	8(%rsi), %rax
addq    %rax, %rax
movq    %rax, 64(%rsp)
mulq 	32(%rsi)
movq 	%rax, %r8
movq 	%rdx, %rbx

movq 	16(%rsi), %rax
addq    %rax, %rax
movq    %rax, 72(%rsp)
mulq 	24(%rsi)
addq 	%rax, %r8
adcq 	%rdx, %rbx

movq 	72(%rsp),%rax
mulq 	32(%rsi)
movq 	%rax, %r10
movq 	%rdx, %r11

movq 	24(%rsi),%rax
mulq 	%rax
addq 	%rax, %r10
adcq 	%rdx, %r11

movq 	24(%rsi),%rax
addq    %rax, %rax
mulq 	32(%rsi)
movq 	%rax, %r12
movq 	%rdx, %r13

movq 	32(%rsi),%rax
mulq 	%rax
movq 	%rax, %r14
movq 	%rdx, %r15

shld    $12, %r8, %rbx
shld    $12, %r10, %r11
shld    $12, %r12, %r13
shld    $12, %r14, %r15

movq    $0xFFFFFFFFFFFFF, %rbp

andq	%rbp, %r8
andq	%rbp, %r10
andq	%rbp, %r12
andq	%rbp, %r14

addq    %rbx, %r10
addq    %r11, %r12
addq    %r13, %r14

movq    $0x1000003D10, %rdi

movq    %r8, %rax
mulq    %rdi
movq    %rax, %r8
movq    %rdx, %r9

movq 	0(%rsi), %rax
mulq 	%rax
addq 	%rax, %r8
adcq 	%rdx, %r9

shld    $12, %r8, %r9

movq    %r10, %rax
mulq    %rdi
movq    %rax, %r10
movq    %rdx, %r11

movq 	0(%rsi), %rax
addq    %rax, %rax
movq    %rax, 72(%rsp)
mulq 	8(%rsi)
addq 	%rax, %r10
adcq 	%rdx, %r11

shld    $12, %r10, %r11

movq    %r12, %rax
mulq    %rdi
movq    %rax, %r12
movq    %rdx, %r13

movq 	72(%rsp),%rax
mulq 	16(%rsi)
addq 	%rax, %r12
adcq 	%rdx, %r13

movq 	8(%rsi),%rax
mulq 	%rax
addq 	%rax, %r12
adcq 	%rdx, %r13

shld    $12, %r12, %r13

movq    %r14, %rax
mulq    %rdi
movq    %rax, %rbx
movq    %rdx, %rbp

movq 	72(%rsp),%rax
mulq 	24(%rsi)
addq 	%rax, %rbx
adcq 	%rdx, %rbp

movq 	64(%rsp),%rax
mulq 	16(%rsi)
addq 	%rax, %rbx
adcq 	%rdx, %rbp

shld    $12, %rbx, %rbp

movq    %r15, %rax
mulq    %rdi
movq    %rax, %r14
movq    %rdx, %r15

movq 	72(%rsp),%rax
mulq 	32(%rsi)
addq 	%rax, %r14
adcq 	%rdx, %r15

movq 	64(%rsp),%rax
mulq 	24(%rsi)
addq 	%rax, %r14
adcq 	%rdx, %r15

movq 	16(%rsi),%rax
mulq 	%rax
addq 	%rax, %r14
adcq 	%rdx, %r15

shld    $16, %r14, %r15

movq    $0xFFFFFFFFFFFFF, %rsi
movq    $0xFFFFFFFFFFFF, %rcx
movq    $0x1000003D1, %rdi

andq	%rsi, %r8
andq	%rsi, %r10
andq	%rsi, %r12
andq	%rsi, %rbx
andq	%rcx, %r14

addq    %r9, %r10
addq    %r11, %r12
addq    %r13, %rbx
addq    %rbp, %r14

movq    %rbx, %rdx
shrq    $52, %rdx
addq    %rdx, %r14
andq	%rsi, %rbx

movq    %r14, %rax
shrq    $48, %rax
addq    %r15, %rax
andq	%rcx, %r14

mulq    %rdi
addq    %rax, %r8
adcq    $0, %rdx
shld    $12, %r8, %rdx
andq	%rsi, %r8
addq    %r10, %rdx

movq    %rdx, %r10
shrq    $52, %rdx
addq    %r12, %rdx
andq	%rsi, %r10

movq    %rdx, %r12
shrq    $52, %rdx
addq    %rbx, %rdx
andq	%rsi, %r12

movq    %rdx, %r13
shrq    $52, %rdx
addq    %rdx, %r14
andq	%rsi, %r13

movq 	56(%rsp), %rdi

movq 	%r8,   0(%rdi)
movq 	%r10,  8(%rdi)
movq 	%r12, 16(%rdi)
movq 	%r13, 24(%rdi)
movq 	%r14, 32(%rdi)

movq 	 8(%rsp), %r12
movq 	16(%rsp), %r13
movq 	24(%rsp), %r14
movq 	32(%rsp), %r15
movq 	40(%rsp), %rbp
movq 	48(%rsp), %rbx
movq 	 0(%rsp), %rsp
ret
.size secp256k1_fe_sqr_inner, .-secp256k1_fe_sqr_inner
