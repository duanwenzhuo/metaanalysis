"""core/exceptions.py — 自定义异常定义。"""

class BioAnalysisError(Exception):
    """Bio-Analysis 系统的基础异常类型"""
    pass


class AnalysisError(BioAnalysisError):
    """分析阶段发生的错误（Phase 4/5）"""
    def __init__(self, message, phase=None, gse_id=None):
        self.phase = phase
        self.gse_id = gse_id
        super().__init__(message)


class DataDownloadError(BioAnalysisError):
    """GEO 数据集下载失败"""
    def __init__(self, gse_id, message):
        self.gse_id = gse_id
        super().__init__(f"[{gse_id}] {message}")


class InvalidHypothesisError(BioAnalysisError):
    """假设文本无效或无法解析"""
    pass


class PhaseError(BioAnalysisError):
    """Phase 执行失败（通用）"""
    def __init__(self, phase, message):
        self.phase = phase
        super().__init__(f"Phase {phase} failed: {message}")


class SkillExecutionError(BioAnalysisError):
    """Skill 执行失败"""
    def __init__(self, skill_name, message):
        self.skill_name = skill_name
        super().__init__(f"Skill '{skill_name}' failed: {message}")


class RegistryError(BioAnalysisError):
    """Registry 读写或状态验证错误"""
    pass


class ConfigurationError(BioAnalysisError):
    """配置缺失或无效（环境变量/配置文件）"""
    pass
