# Control code for the paper "Aggregate Grid-Forming Control of Heterogeneous Distributed Energy Resources"
# Optimization code for the paper "Optimal Dynamic Ancillary Services Provision Based on Local Power Grid Perception"

This repository contains the simulation files to reproduce the results of the paper https://arxiv.org/abs/2410.14912; it also contains the optimization code files to reproduce the results of the paper https://arxiv.org/abs/2401.17793

**MATLAB Version: 2022b**

**He, Xiuqiang, Josué Duarte, Verena Häberle, and Florian Dörfler. "Aggregate Grid-Forming Control of Heterogeneous Distributed Energy Resources." arXiv preprint arXiv:2410.14912 (2024).**

**Abstract:** This article presents an aggregate grid-forming control for heterogeneous distributed energy resources (DERs). The proposed control achieves a desired aggregate grid-forming response by coordinating power contributions among multiple different DERs. Unlike existing aggregate control strategies that are typically objective-specific or topology-specific, this article proposes a generic, flexible, and modular control design. The design supports four basic module types---AC- or DC-coupling and AC- or DC-output topological arrangements---adequately accommodating diverse DER integration scenarios, e.g., AC, DC, AC/DC hybrid microgrids, hybrid energy storage systems, or hybrid renewable power plants. The grid-forming control design is systematically developed by aggregating DER dynamics and disaggregating control objectives for the four basic modules, then extended to modular configurations through a top-down approach. The grid-forming performance is comprehensively validated through simulation. This modular control design provides scalable and standardizable grid interfaces, enabling effective aggregation of heterogeneous DERs to jointly achieve grid-forming control and operations.

**Verena Häberle, Xiuqiang He, Linbin Huang, Eduardo Prieto-Araujo, Florian Dörfler. "Optimal Dynamic Ancillary Services Provision Based on Local Power Grid Perception." arXiv preprint arXiv:2401.17793 (2024).**

**Abstract:** In this paper, we propose a systematic closed-loop approach to provide optimal dynamic ancillary services with converter-interfaced generation systems based on local power grid perception. In particular, we structurally encode dynamic ancillary services such as fast frequency and voltage regulation in the form of a parametric transfer function matrix, which includes several parameters to define a set of different feasible response behaviors, among which we aim to find the optimal one to be realized by the converter system. Our approach is based on a so-called "perceive-and-optimize" (P&O) strategy: First, we identify a grid dynamic equivalent at the interconnection terminals of the converter system. Second, we consider the closed-loop interconnection of the identified grid equivalent and the parametric transfer function matrix, which we optimize for the set of transfer function parameters, resulting in a stable and optimal closed-loop performance for ancillary services provision. In the process, we ensure that grid-code and device-level requirements are satisfied. Finally, we demonstrate the effectiveness of our approach in different numerical case studies based on a modified Kundur two-area test system.


**Bibtex to cite the papers**

@article{he2024aggregate,
  title={Aggregate Grid-Forming Control of Heterogeneous Distributed Energy Resources},
  author={He, Xiuqiang and Duarte, Josu{\'e} and H{\"a}berle, Verena and D{\"o}rfler, Florian},
  journal={arXiv preprint arXiv:2410.14912},
  year={2024}
}

@article{häberle2024optimaldynamicancillaryservices,
      title={Optimal Dynamic Ancillary Services Provision Based on Local Power Grid Perception}, 
      author={Verena Häberle and Xiuqiang He and Linbin Huang and Eduardo Prieto-Araujo and Florian Dörfler},
      year={2024},
      eprint={2401.17793},
      archivePrefix={arXiv},
      primaryClass={eess.SY},
      url={https://arxiv.org/abs/2401.17793}, 
}