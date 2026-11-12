# D-STAR (AMBE 2400) reference harness, record format in test/dstar_ref.h
set(XOMBE_DSTAR_REFS ${CMAKE_CURRENT_SOURCE_DIR}/vectors/dstar)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dstar_out)  # outputs, also of the tool tests that need no Python
find_package(Python3 COMPONENTS Interpreter QUIET)

add_library(xombe_dstar_ref STATIC test/dstar_ref.c)
xombe_target_options(xombe_dstar_ref)
xombe_executable(test_dstar_refs test/test_dstar_refs.c)
target_link_libraries(test_dstar_refs xombe_dstar_ref)
add_test(NAME dstar_refs COMMAND test_dstar_refs ${XOMBE_DSTAR_REFS})

# xombe_dstar_test(<name> [DIR <set>] [CLIPS <clip>...] INPUT <file-or-suffix> REF <file-or-suffix>
#                  RUN <program> [<args>...] [CMP <dstar_cmp.py options>...])
#
# Adds one ctest per clip (or one test without CLIPS), named dstar_<name>[_<clip>], that runs
#   <program> <args> <set>/<clip><INPUT> <out>
# and then compares <out> with <set>/<clip><REF> using test/dstar_cmp.py <CMP options>; a REF ending
# in .raw is compared as PCM (--pcm).  <set> is vectors/dstar unless DIR gives another.  <program> may be a
# target name.  The program reads the input (a stream, or a reference file it takes its inputs from) and writes
# its dump with test/dstar_ref.h's writer, using the reference file's column names.  Needs Python 3 (standard
# library only): without it each test is a placeholder that ctest reports as skipped (label no_python).
function(xombe_dstar_test name)
  cmake_parse_arguments(A "" "DIR;INPUT;REF" "CLIPS;RUN;CMP" ${ARGN})
  if(A_DIR)
    set(dir ${A_DIR})
  else()
    set(dir ${XOMBE_DSTAR_REFS})
  endif()
  list(GET A_RUN 0 prog)
  list(REMOVE_AT A_RUN 0)
  if(TARGET ${prog})
    set(prog $<TARGET_FILE:${prog}>)
  endif()
  set(pcm "")
  if(A_REF MATCHES "\\.raw$")
    set(pcm --pcm)
  endif()
  if(NOT A_CLIPS)
    set(A_CLIPS "-")
  endif()
  foreach(clip ${A_CLIPS})
    if(clip STREQUAL "-")
      set(clip "")
      set(tname dstar_${name})
    else()
      set(tname dstar_${name}_${clip})
    endif()
    if(NOT Python3_Interpreter_FOUND)
      add_test(NAME ${tname} COMMAND sh -c "exit 77")
      set_tests_properties(${tname} PROPERTIES SKIP_RETURN_CODE 77 LABELS no_python)
      set_property(GLOBAL APPEND PROPERTY XOMBE_SKIPPED_TESTS ${tname})
      continue()
    endif()
    set(out ${CMAKE_CURRENT_BINARY_DIR}/dstar_out/${tname}.out)
    add_test(NAME ${tname} COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/dstar_cmp.py
             ${pcm} ${A_CMP} ${dir}/${clip}${A_REF} ${out} --run ${prog} ${A_RUN} ${dir}/${clip}${A_INPUT} ${out})
  endforeach()
endfunction()

# In-repo clips: sp (speech), nz (speech in noise), tn (tones), kx (Knox tones, near silence) and
# sper (sp with special frames and bit errors; D-Star input only); rate-33 inputs: sp, tn, kx, sper, tx.
set(XOMBE_DSTAR_CLIPS sp nz tn kx sper)
set(XOMBE_DSTAR_R33_CLIPS sp tn kx sper tx)

# The comparison machinery itself: each reference against a copy of itself.
xombe_dstar_test(selfcheck CLIPS sp sper INPUT .fec.txt.gz REF .fec.txt.gz RUN ${CMAKE_COMMAND} -E copy)
xombe_dstar_test(selfcheck_pcm CLIPS sper INPUT .dec.raw REF .dec.raw RUN ${CMAKE_COMMAND} -E copy)
xombe_dstar_test(selfcheck_quant INPUT quant_closed.txt.gz REF quant_closed.txt.gz RUN ${CMAKE_COMMAND} -E copy
                 CMP --key tset,parent,pass)

# deq floats: float32 against the model's float64, measured max |error|:
# a 3.2e-6 (log2; float32 ulp is 9.5e-7 for |a| in [8, 16)), Tv 1.0e-6, G 4.4e-7, f0 1.8e-9.
set(XOMBE_DSTAR_DEQ_TOL --tol *=1e-5:0)
if(XOMBE_DSTAR)
  # P1: FEC, fields, the error rule and the encoder's packer (exact), and dequantization against the mirror.
  xombe_internal_executable(test_dstar_tables test/test_dstar_tables.c)
  target_link_libraries(test_dstar_tables xombe_dstar_ref)
  add_test(NAME dstar_tables COMMAND test_dstar_tables ${XOMBE_DSTAR_REFS})
  xombe_internal_executable(dstar_dump_fec test/dstar_dump_fec.c)
  target_link_libraries(dstar_dump_fec xombe_dstar_ref)
  xombe_internal_executable(dstar_dump_deq test/dstar_dump_deq.c)
  target_link_libraries(dstar_dump_deq xombe_dstar_ref)
  xombe_dstar_test(fec CLIPS ${XOMBE_DSTAR_CLIPS} INPUT .dstar REF .fec.txt.gz RUN dstar_dump_fec)
  xombe_dstar_test(fec_enc INPUT fec_enc.txt.gz REF fec_enc.txt.gz RUN dstar_dump_fec --enc
                   CMP --tol level=0:0)
  xombe_dstar_test(deq CLIPS ${XOMBE_DSTAR_CLIPS} INPUT .dstar REF .deq.txt.gz RUN dstar_dump_deq
                   CMP ${XOMBE_DSTAR_DEQ_TOL})
  # P2: the decoder (dstar_dump_dec).
  xombe_internal_executable(dstar_dump_dec test/dstar_dump_dec.c)
  target_link_libraries(dstar_dump_dec xombe_dstar_ref)
  # tolerances, measured:
  # dec: float32 against the model's float64: Pia 7.2e-6, PiG 3.5e-6, lam 4.6e-6 (log2 units), Pif0 3.9e-7
  #   (the initial memory's f0); fbar 0.0033 Hz (AMBE+2's initial average pitch, 111.63 Hz against 111.6333).
  # half: a 5.2e-5 at |a| ~ 12 (6 printed digits), M 3e-5 relative (2^a), rms 0.5 + 0.2%.
  # PCM: the AMBE+2 synthesis' own C-vs-Python difference (relative ~1e-3, max 5 LSB on AMBE+2's es2). tn and sper
  #   meet the gate up to their first dual tone (frames 25 and 31), then drift in phase: xombe's tone table f0
  #   (US 7,970,606 B2) and the model's fitted one differ by up to 0.007 Hz. pcm_all bounds the whole clip
  set(XOMBE_DSTAR_DEC_TOL --tol *=1e-5:0 --tol fbar=0.01:0)
  set(XOMBE_DSTAR_HALF_TOL --tol *=1e-5:1e-5 --tol M=1e-5:3e-5 --tol rms=0.5:2e-3)
  set(XOMBE_DSTAR_PCM_GATE --lsb 4 --max-rms 3 --min-within 0.995)
  xombe_dstar_test(dec CLIPS ${XOMBE_DSTAR_CLIPS} INPUT .dstar REF .dec.txt.gz RUN dstar_dump_dec
                   CMP ${XOMBE_DSTAR_DEC_TOL})
  xombe_dstar_test(half CLIPS sp nz tn kx sper INPUT .dstar REF .half.txt.gz RUN dstar_dump_dec --half
                   CMP ${XOMBE_DSTAR_HALF_TOL})
  xombe_dstar_test(pcm CLIPS sp nz kx INPUT .dstar REF .dec.raw RUN dstar_dump_dec --pcm
                   CMP ${XOMBE_DSTAR_PCM_GATE})
  xombe_dstar_test(pcm CLIPS tn INPUT .dstar REF .dec.raw RUN dstar_dump_dec --pcm
                   CMP ${XOMBE_DSTAR_PCM_GATE} --samples 0:4000)
  xombe_dstar_test(pcm CLIPS sper INPUT .dstar REF .dec.raw RUN dstar_dump_dec --pcm
                   CMP ${XOMBE_DSTAR_PCM_GATE} --samples 0:4960)
  xombe_dstar_test(pcm_all CLIPS tn INPUT .dstar REF .dec.raw RUN dstar_dump_dec --pcm
                   CMP --lsb 4 --max-rms 110 --min-within 0.55)
  xombe_dstar_test(pcm_all CLIPS sper INPUT .dstar REF .dec.raw RUN dstar_dump_dec --pcm
                   CMP --lsb 4 --max-rms 4 --min-within 0.94)
endif()

# P3: the closed-loop amplitude quantizer (src/common/cquant.c) on the quantizer cases.
if(XOMBE_DSTAR AND XOMBE_ENCODER)
  xombe_internal_executable(dstar_dump_quant test/dstar_dump_quant.c)
  target_link_libraries(dstar_dump_quant xombe_dstar_ref)
  set(XOMBE_DSTAR_QUANT_CMP --tol *=1e-5:1e-4)
  xombe_dstar_test(quant_closed INPUT quant_closed.txt.gz REF quant_closed.txt.gz RUN dstar_dump_quant closed
                   CMP --key tset,parent,pass ${XOMBE_DSTAR_QUANT_CMP})
  xombe_dstar_test(quant_dstar INPUT quant_dstar.txt.gz REF quant_dstar.txt.gz RUN dstar_dump_quant dstar
                   CMP --key id ${XOMBE_DSTAR_QUANT_CMP})
  xombe_dstar_test(quant_ambe2 INPUT quant_ambe2.txt.gz REF quant_ambe2.txt.gz RUN dstar_dump_quant ambe2
                   CMP --key id ${XOMBE_DSTAR_QUANT_CMP})
  # capped at 100 units (the encoders use XOMBE_CQ_CAP): the bound per search, the summed error within 1e-3
  # of the exact search's (checked by dstar_dump_quant), and 99% of the calls as the reference
  set(cmin_closed --min b=0.99 --min err=0.99)
  set(cmin_dstar --min b=0.99 --min ahat=0.99 --min err=0.99)
  set(cmin_ambe2 --min bclosed=0.99 --min pick=0.99 --min b=0.99 --min ahat=0.99 --min err=0.99)
  set(ckey_closed tset,parent,pass)
  set(ckey_dstar id)
  set(ckey_ambe2 id)
  foreach(k closed dstar ambe2)
    xombe_dstar_test(quant_${k}_cap INPUT quant_${k}.txt.gz REF quant_${k}.txt.gz RUN dstar_dump_quant --cap 100 ${k}
                     CMP --key ${ckey_${k}} ${XOMBE_DSTAR_QUANT_CMP} ${cmin_${k}})
  endforeach()
endif()

# P4b: the encoder back end (src/dstar/dstar_enc.c) on the model's recorded inputs (test/dstar_dump_enc_be.c).
# Exact: b0 from the refined period, b1 from the recorded a and w, the fields and bytes of the recorded targets
# quantized against the mirror over the stream. m_noise is the generated table (16 bits, relative error < 3e-5).
# --spec: the voicing on the C spectrum of the input (C high-pass, float), against the model's.
if(XOMBE_DSTAR AND XOMBE_ENCODER)
  xombe_internal_executable(dstar_dump_enc_be test/dstar_dump_enc_be.c)
  target_link_libraries(dstar_dump_enc_be xombe_dstar_ref)
  set(XOMBE_DSTAR_ENC_CLIPS sp nz tn kx)
  xombe_dstar_test(enc_be CLIPS ${XOMBE_DSTAR_ENC_CLIPS} INPUT .enc.txt.gz REF .enc.txt.gz RUN dstar_dump_enc_be
                   CMP --tol mn=0:1e-4)
  # measured: a within 0.01 on 99.84-100% of frames (worst 0.073), b1 equal on all but one frame, mv and w differ
  # by the input high-pass (4-digit coefficients against scipy's Butterworth)
  set(XOMBE_DSTAR_ENC_SPEC_CMP --tol mv=0:0.01 --tol w=0:0.01 --tol a=0.01:0 --min mv=0.9 --min w=0.95 --min a=0.995
      --min b1=0.999)
  # The q15 kernels' block-floating-point FFT (~60-65 dB) moves the metric further: measured minimums over the
  # same clips w 76.7%, mv 79.2%, a 94.8%, b1 99.77%
  if(XOMBE_KERNELS STREQUAL "q15")
    set(XOMBE_DSTAR_ENC_SPEC_CMP --tol mv=0:0.01 --tol w=0:0.01 --tol a=0.01:0 --min mv=0.78 --min w=0.75
        --min a=0.94 --min b1=0.997)
  endif()
  xombe_dstar_test(enc_be_spec CLIPS ${XOMBE_DSTAR_ENC_CLIPS} INPUT .enc.txt.gz REF .enc.txt.gz
                   RUN dstar_dump_enc_be --spec CMP ${XOMBE_DSTAR_ENC_SPEC_CMP})

  # P4c: the whole encoder with the D-Star output on the clip's PCM (test/dstar_dump_enc.c): the tone pipeline's
  # detection, run rule and final tone, and the frame class, exact against the model; levels within 0.02.
  xombe_internal_executable(dstar_dump_enc test/dstar_dump_enc.c)
  target_link_libraries(dstar_dump_enc xombe_dstar_ref)
  set(XOMBE_DSTAR_ENC_TONE_CMP --tol det=0.02:0 --tol trk=0.02:0 --tol tone=0.02:0 --skip b0,b1,out)
  # kx's Knox frames take their pairs from xombe's tone table (US 7,970,606 B2's f0 and harmonics); the model's
  # are fitted, and every kx frame still matches, bytes included.
  xombe_dstar_test(enc CLIPS sp nz tn kx INPUT .in.raw REF .enc.txt.gz RUN dstar_dump_enc CMP ${XOMBE_DSTAR_ENC_TONE_CMP})
  xombe_dstar_test(enc_dual CLIPS tn INPUT .in.raw REF .enc.txt.gz RUN dstar_dump_enc --dual
                   CMP ${XOMBE_DSTAR_ENC_TONE_CMP})
endif()

# P5: the rate-33 -> D-STAR transcoder (src/dstar/xcode.c) on the rate-33 streams (test/dstar_dump_xc.c).
# Floats: f0, a, G, soft and tgt within 1e-5 (float32 against the model's float64). Knox frames (kx, tx) take f0
# and harmonics from xombe's tone table (US 7,970,606 B2), the model its fitted one: b0 differs for IDs 145, 150,
# 151, 153, 154 and 157, and the mirror carries the difference on
if(XOMBE_DSTAR)
  xombe_internal_executable(dstar_dump_xc test/dstar_dump_xc.c)
  target_link_libraries(dstar_dump_xc xombe_dstar_ref)
  set(XOMBE_DSTAR_XC_CMP --tol *=1e-5:1e-5)
  set(XOMBE_DSTAR_XC_KNOX --min tgt=0.4 --min tb0=0.44 --min b=0.17 --min out=0.17)
  xombe_dstar_test(xc CLIPS sp tn sper INPUT .r33 REF .xc.txt.gz RUN dstar_dump_xc CMP ${XOMBE_DSTAR_XC_CMP})
  xombe_dstar_test(xc CLIPS kx tx INPUT .r33 REF .xc.txt.gz RUN dstar_dump_xc CMP ${XOMBE_DSTAR_XC_CMP}
                   ${XOMBE_DSTAR_XC_KNOX})
  # the tool's stream equals the model's output
  if(XOMBE_BUILD_TOOLS)
    foreach(clip sp tn sper)
      add_test(NAME dstar_xc_tool_${clip} COMMAND sh -c "\"$1\" \"$2\" \"$3\" && \"$4\" -E compare_files \"$3\" \"$5\""
               xc $<TARGET_FILE:xombe-xc> ${XOMBE_DSTAR_REFS}/${clip}.r33
               ${CMAKE_CURRENT_BINARY_DIR}/dstar_out/xc_tool_${clip}.dstar ${CMAKE_COMMAND}
               ${XOMBE_DSTAR_REFS}/${clip}.xc.dstar)
    endforeach()
  endif()
endif()

# the D-STAR -> rate-33 transcoder (src/dstar/xrev.c) on the D-STAR streams (test/dstar_dump_xrev.c).
# floats: f0, a and tgt within 1e-5 (float32 against the model's float64)
if(XOMBE_DSTAR AND XOMBE_ENCODER)
  xombe_internal_executable(dstar_dump_xrev test/dstar_dump_xrev.c)
  target_link_libraries(dstar_dump_xrev xombe_dstar_ref)
  set(XOMBE_DSTAR_XREV_CMP --tol *=1e-5:1e-5)
  xombe_dstar_test(xrev CLIPS sp tn kx sper INPUT .dstar REF .xrev.txt.gz RUN dstar_dump_xrev
                   CMP ${XOMBE_DSTAR_XREV_CMP})
  # with the transcoder's XOMBE_CQ_CAP (the tool's search): the same records
  xombe_dstar_test(xrev_capped CLIPS sp tn kx sper INPUT .dstar REF .xrev.txt.gz RUN dstar_dump_xrev --capped
                   CMP ${XOMBE_DSTAR_XREV_CMP})
  if(XOMBE_BUILD_TOOLS)
    foreach(clip sp tn kx sper)
      add_test(NAME dstar_xrev_tool_${clip}
               COMMAND sh -c "\"$1\" --reverse \"$2\" \"$3\" && \"$4\" -E compare_files \"$3\" \"$5\""
               xrev $<TARGET_FILE:xombe-xc> ${XOMBE_DSTAR_REFS}/${clip}.dstar
               ${CMAKE_CURRENT_BINARY_DIR}/dstar_out/xrev_tool_${clip}.r33 ${CMAKE_COMMAND}
               ${XOMBE_DSTAR_REFS}/${clip}.xrev.r33)
    endforeach()
  endif()
endif()
