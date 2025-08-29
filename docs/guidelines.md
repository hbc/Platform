# Data Storage, Pipelines, and Analysis Environment Guidelines

**Audience:** Internal team (Shannan, Lorena, Alex, Noor, Elizabeth, Ruitong, Zhu, James, Upen, Will, Emma, others)
**Purpose:** Establish consistent practices for data storage, pipeline execution, and downstream analysis across FAS, O2, MIX, AWS, and local environments.

---

## 1. Data Storage

### 1.1 Storage Locations

* **FAS**

  * For FAS/HSPH-affiliated projects.
  * Intended use: *final* and *downstream objects*.
  * Scratch storage: temporary data only; May be auto-synced to O2 for non-FAS projects.

* **O2**

  * For HMS-affiliated and non FAS/HSPH projects.

* **Globus**

  * Preferred for FAS–O2 syncing (benchmark: 1 TB in \~15 minutes).

### 1.2 Data Lifecycle Policy

* **Raw Data**:

  * Work only in scratch
  * Used only for pipelines
  * Returned to clients after processing
  * Explore “cold storage” with retrieval fees

* **Pipeline Outputs Data**:

  * Long-term retention on FAS.
  * Accessible to clients/collaborators.

* **Downstream Objects**:

  * Retained on FAS for reproducibility and reanalysis.

### 1.3 Data Management and Data Flow Practices

* Maintain **strict folder structure** for consistency. Follow project names in **Trello Cards**
* Ensures all PI folder directories have `group r+w` permissions.
* Every analyst works in **their workspace**:
    - It could be FAS scratch or FAS user space
    - **First step is to clone repositiory** (repository would be ready for analysts to clone)
    - Data (Primary-pipeline outputs, Secondary-files and objects) **always** in project folder . 
        - Easy level: use full paths to project folder
        - Advance level: use symlinks if it is easy for you. Add that step into the readme if you do use symlinks.

This ensure multiple people working in the same project, avoid GitHub issues, ensure good data management practices.

## 2. Computational Environments

### 2.1 Pipelines (Primary Analysis)

* **FAS**: Default environment for most projects.
    * Monitor with Seqera (Alex to assist).
* **Seqera + AWS**: For small workloads (< 20 samples) or projects requiring cloud scalability.

### 2.2 Downstream Analysis (Secondary / Exploratory)

* **Local**: Preferred by many researchers for lightweight analysis.
    * Downstream objects need to be put back in project directory at FAS or O2
* **FAS**: Standard for large datasets, reproducibility, and shared work.
* **O2**: Used for HMS projects and where performance is sufficient.

## 3. Open Questions / Action Items

1. **Storage Cost Policy**

   * Should clients pay for raw data retention?
   * Should we implement cold storage + retrieval fees?

2. **Permissions and Automation**

   * Confirm with O2 team about cron jobs for PI folders (`group r+w`).
   * Define rules for excluding specific directories in sync.

3. **Team Coordination**

   * Shannan & Lorena: schedule meeting with FASRC to discuss quotas, automation, and permissions.
